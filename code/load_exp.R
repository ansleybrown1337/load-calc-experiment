
# load-calc-experiment
# Script to calculate water analyte loads and showcase limitations of methods used
# Updated:
#   - robust ISCO flow import (variable header length)
#   - optional volume estimation from avg_flow_gpm when volume_gal is missing/NA
#   - Plotly helper for interactive flow review + bottle-colored sample-event markers
#   - Fe/Se unit normalization: WQ for Fe/Se assumed ug/L, converted to mg/L before load calc
#   - TN derived analyte (TKN + NO3_N + NO2_N)
#   - Uncertainty handling: if SD cannot be estimated (n<2), bounds are NA (not 0)
#
# IMPORTANT: All loads are returned in kg (for all analytes, including Fe/Se).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

# ============================ ANALYTE ABBREVIATION DICTIONARY ============================

analyteAbbr.dict <- list(
  "TKN"   = c("Nitrogen, Total Kjeldahl"),
  "NO2_N" = c("Nitrogen, Nitrite  (As N)", "NITRITE AS N"),
  "PO4_P" = c("Phosphorus, Total Orthophosphate (as P)", "ORTHOPHOSPHATE AS P"),
  "TP"    = c("Phosphorus, Total (As P)", "TOTAL PHOSPHORUS"),
  "TDS"   = c("Total Dissolved Solids (Residue, Filterable)"),
  "NO3_N" = c("Nitrogen, Nitrate (As N)", "NITRATE AS N"),
  "TSS"   = c("Suspended Solids (Residue, Non-Filterable)", "TSS"),
  "Fe"    = c("Iron, Total"),
  "Se"    = c("Selenium, Total", "SELENIUM"),
  "pH"    = c("pH", "POTENTIAL HYDROGEN"),
  "EC25"  = c("Specific Conductance", "ELECTRICAL CONDUCTIVITY")
)

legacy_analyte_label_map <- c(
  "Fe (g not kg)" = "Fe",
  "Se (g not kg)" = "Se"
)

# ============================ SMALL HELPERS ============================

clean_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))
}

# NOTE: ND and values like "<0.01" are treated as ZERO (per your preference)
parse_result_to_numeric_nd_as_zero <- function(x) {
  x_chr <- trimws(as.character(x))
  is_lt <- grepl("^<", x_chr)
  nums  <- suppressWarnings(as.numeric(x_chr))
  out   <- ifelse(is_lt, 0, nums)      # <MDL -> 0
  out[grepl("^\\s*(ND|nd|non[- ]?detect(ed)?|n/?d)\\s*$", x_chr)] <- 0
  out[!is.finite(out)] <- NA_real_
  out
}

normalize_analyte_abbr <- function(x) {
  x_chr <- as.character(x)
  x_chr <- ifelse(x_chr %in% names(legacy_analyte_label_map), legacy_analyte_label_map[x_chr], x_chr)
  x_chr
}

# Fe and Se are assumed reported in ug/L in the WQ file; convert ug/L -> mg/L
normalize_wq_units_to_mgL <- function(result, analyte_abbr) {
  out <- result
  is_fe_se <- analyte_abbr %in% c("Fe", "Se")
  out[is_fe_se & is.finite(out)] <- out[is_fe_se & is.finite(out)] / 1000
  out
}

# ============================ IMPORT DATA ============================

import_wq_data <- function(file_path) {
  read.csv(file_path, stringsAsFactors = TRUE)
}

# Robust ISCO flow import: auto-detect first datetime row, then return table portion
import_flow_data_flex <- function(file_path, tz = "America/Denver") {
  raw <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (ncol(raw) < 2) stop("Flow file appears to have <2 columns. Check delimiter / file.", call. = FALSE)
  
  dt_try <- suppressWarnings(lubridate::parse_date_time(
    raw[[1]],
    orders = c("mdY HMS p", "mdY HM p", "mdY HMS", "mdY HM"),
    tz = tz
  ))
  
  first_data_row <- which(!is.na(dt_try))[1]
  
  if (is.na(first_data_row)) {
    stop(
      paste0(
        "Could not find any rows whose first column parses as a datetime. ",
        "Expected something like 'm/d/Y h:m:s AM/PM' or 'm/d/Y H:M[:S]'."
      ),
      call. = FALSE
    )
  }
  
  raw[first_data_row:nrow(raw), , drop = FALSE]
}
# Backward-compatible alias
import_flow_data <- function(file_path) {
  import_flow_data_flex(file_path, tz = "America/Denver")
}

# ============================ PROCESS DATA ============================

process_wq_data <- function(wq_data, treatment_filter) {
  wq_data %>%
    mutate(
      collected = as.POSIXct(parse_date_time(
        collected, orders = c("mdY HM","mdY","Ymd HMS","Ymd","d b Y HM")
      ), tz = "America/Denver"),
      received = as.POSIXct(parse_date_time(
        received, orders = c("mdY HM","mdY","Ymd HMS","Ymd","d b Y HM")
      ), tz = "America/Denver"),

      treatment.name = as.character(treatment.name),

      analyte_abbr = vapply(analyte, function(a) {
        match <- names(analyteAbbr.dict)[sapply(analyteAbbr.dict, function(x) a %in% x)]
        if (length(match) == 0) NA_character_ else match
      }, FUN.VALUE = character(1)),

      result = parse_result_to_numeric_nd_as_zero(result),

      analyte_abbr = normalize_analyte_abbr(analyte_abbr),
      result = normalize_wq_units_to_mgL(result, analyte_abbr)
    ) %>%
    {
      if (all(is.na(.$treatment.name))) {
        stop("All values in 'treatment.name' are NA. Verify your dataset.", call. = FALSE)
      } else {
        filter(., treatment.name == treatment_filter | event.type == treatment_filter)
      }
    } %>%
    {
      if (any(is.na(.$analyte_abbr))) {
        stop("Some analytes could not be mapped. Update analyteAbbr.dict accordingly.", call. = FALSE)
      } else {
        .
      }
    }
}

process_flow_data_flex <- function(flow_data, tz = "America/Denver") {
  num_columns <- ncol(flow_data)
  
  if (num_columns == 8) {
    column_names <- c("datetime","min_flow_raw","min_time","max_flow_gpm",
                      "max_time","avg_flow_gpm","volume_gal","sample_event")
  } else if (num_columns == 6) {
    column_names <- c("datetime","min_flow_raw","max_flow_gpm",
                      "avg_flow_gpm","volume_gal","sample_event")
  } else {
    stop(
      paste0(
        "Unexpected number of columns in flow data: ", num_columns, ". ",
        "You may have passed the wrong file, or the export format differs."
      ),
      call. = FALSE
    )
  }
  
  out <- flow_data %>%
    setNames(column_names) %>%
    mutate(
      datetime = suppressWarnings(lubridate::parse_date_time(
        datetime,
        orders = c("mdY HMS p", "mdY HM p", "mdY HMS", "mdY HM"),
        tz = tz
      )),
      across(
        intersect(c("min_flow_raw","max_flow_gpm","avg_flow_gpm","volume_gal","sample_event"), names(.)),
        clean_numeric
      )
    ) %>%
    # Critical: enforce time order to prevent Plotly “looping”
    arrange(datetime)
  
  if ("min_time" %in% names(out)) {
    out <- out %>%
      mutate(min_time = as.POSIXct(
        paste(as.Date(datetime), min_time),
        format = "%Y-%m-%d %I:%M:%S %p", tz = tz
      ))
  }
  if ("max_time" %in% names(out)) {
    out <- out %>%
      mutate(max_time = as.POSIXct(
        paste(as.Date(datetime), max_time),
        format = "%Y-%m-%d %I:%M:%S %p", tz = tz
      ))
  }
  
  if (!("min_flow_gpm" %in% names(out))) {
    out <- out %>% mutate(min_flow_gpm = if ("min_flow_raw" %in% names(out)) min_flow_raw else NA_real_)
  }
  
  # (keep your existing volume estimation logic here unchanged, if present)
  need_est <- (!("volume_gal" %in% names(out))) || all(is.na(out$volume_gal))
  if (need_est) {
    ok <- is.finite(out$avg_flow_gpm) & !is.na(out$datetime)
    dt_ok <- out$datetime[ok]
    
    if (length(dt_ok) >= 2) {
      delta_mins <- as.numeric(difftime(dplyr::lead(out$datetime), out$datetime, units = "mins"))
      
      diffs <- as.numeric(difftime(dt_ok[-1], dt_ok[-length(dt_ok)], units = "mins"))
      diffs <- diffs[is.finite(diffs) & diffs > 0]
      med_step <- if (length(diffs) > 0) stats::median(diffs) else 15
      
      delta_mins[!is.finite(delta_mins) | delta_mins <= 0] <- med_step
      out$volume_gal <- out$avg_flow_gpm * delta_mins
    } else {
      out$volume_gal <- NA_real_
    }
  }
  
  out %>% mutate(volume_L = ifelse(is.finite(volume_gal), volume_gal * 3.78541, NA_real_))
}

process_flow_data <- function(flow_data) {
  process_flow_data_flex(flow_data, tz = "America/Denver")
}

# ============================ AGGREGATE FLOW DATA ============================

aggregate_flow_data <- function(flow_data, user_interval = 4) {
  flow_data %>%
    mutate(
      time_block = as.POSIXct(
        floor(as.numeric(datetime) / (user_interval * 3600)) * (user_interval * 3600),
        origin = "1970-01-01", tz = "America/Denver"
      )
    ) %>%
    group_by(time_block) %>%
    summarise(
      avg_flow_gpm = mean(avg_flow_gpm, na.rm = TRUE),
      volume_L = if ("volume_L" %in% names(flow_data) && any(!is.na(flow_data$volume_L))) {
        sum(volume_L, na.rm = TRUE)
      } else {
        NA_real_
      },
      sample_event = if ("sample_event" %in% names(flow_data)) sum(sample_event, na.rm = TRUE) else NA_real_,
      min_flow_gpm = if ("min_flow_gpm" %in% names(flow_data) && any(!is.na(flow_data$min_flow_gpm))) {
        suppressWarnings(min(min_flow_gpm, na.rm = TRUE))
      } else {
        NA_real_
      },
      max_flow_gpm = if ("max_flow_gpm" %in% names(flow_data) && any(!is.na(flow_data$max_flow_gpm))) {
        suppressWarnings(max(max_flow_gpm, na.rm = TRUE))
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      min_flow_gpm = ifelse(is.infinite(min_flow_gpm), NA_real_, min_flow_gpm),
      max_flow_gpm = ifelse(is.infinite(max_flow_gpm), NA_real_, max_flow_gpm)
    ) %>%
    select(where(~ !all(is.na(.))))
}

# ============================ COMPUTE ANALYTE SUMMARY ============================

# If SD cannot be estimated (fewer than 2 finite results), sd is NA (not 0)
compute_analyte_summary <- function(wq_data) {
  wq_data %>%
    group_by(analyte_abbr) %>%
    summarise(
      n_finite  = sum(is.finite(result)),
      avg_result = mean(result, na.rm = TRUE),
      sd_result  = ifelse(n_finite >= 2, sd(result, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    split(.$analyte_abbr) %>%
    lapply(function(df) list(avg = df$avg_result, sd = df$sd_result))
}

# Add TN as a derived analyte in the *summary* space (mg/L)
add_total_nitrogen_to_summary <- function(analyte_summary,
                                          tn_name = "TN",
                                          components = c("TKN", "NO3_N", "NO2_N"),
                                          drop_components = FALSE) {

  if (!all(components %in% names(analyte_summary))) {
    missing <- components[!(components %in% names(analyte_summary))]
    message("TN not created, missing components: ", paste(missing, collapse = ", "))
    return(analyte_summary)
  }

  mu <- sum(vapply(components, function(k) analyte_summary[[k]]$avg, numeric(1)), na.rm = TRUE)

  sds <- vapply(components, function(k) analyte_summary[[k]]$sd, numeric(1))
  sd_ok <- all(is.finite(sds))
  sd_tn <- if (sd_ok) sqrt(sum(sds^2)) else NA_real_

  analyte_summary[[tn_name]] <- list(avg = mu, sd = sd_tn)

  if (isTRUE(drop_components)) {
    analyte_summary[components] <- NULL
  }

  analyte_summary
}

# ============================ LOAD CALCULATION ============================

compute_load <- function(analyte_obj, analyte_name, flow_data) {
  if (!analyte_name %in% names(analyte_obj)) return(NA)

  avg_conc     <- analyte_obj[[analyte_name]]$avg   # mg/L
  sd_conc      <- analyte_obj[[analyte_name]]$sd    # mg/L (may be NA)
  total_volume <- sum(flow_data$volume_L, na.rm = TRUE)  # L

  if (!is.finite(total_volume) || total_volume <= 0) {
    warning("Total flow volume is non-finite or <= 0; loads returned as NA for analyte: ", analyte_name, call. = FALSE)
    return(data.frame(
      analyte = analyte_name,
      load_kg = NA_real_, load_upper_kg = NA_real_, load_lower_kg = NA_real_
    ))
  }

  # Point estimate always computed if avg_conc is finite
  load_kg <- avg_conc * total_volume / 1e6

  # If SD is not available, keep bounds as NA (do NOT show 0)
  if (!is.finite(sd_conc)) {
    return(data.frame(
      analyte = analyte_name,
      load_kg = load_kg,
      load_upper_kg = NA_real_,
      load_lower_kg = NA_real_
    ))
  }

  upper <- (avg_conc + sd_conc)
  lower <- pmax(0, avg_conc - sd_conc)

  data.frame(
    analyte        = analyte_name,
    load_kg        = load_kg,
    load_upper_kg  = upper    * total_volume / 1e6,
    load_lower_kg  = lower    * total_volume / 1e6
  )
}

compute_all_loads <- function(analyte_obj, flow_data, drop_analytes = c("EC25", "pH")) {
  analyte_names <- setdiff(names(analyte_obj), drop_analytes)

  if (length(drop_analytes) > 0) {
    message("Note: The following analytes were excluded from load calculations because they are not mass-based: ",
            paste(drop_analytes, collapse = ", "))
  }

  if (length(analyte_names) == 0) {
    return(dplyr::tibble(
      analyte = character(),
      load_kg = double(),
      load_upper_kg = double(),
      load_lower_kg = double()
    ))
  }

  loads <- lapply(analyte_names, function(analyte) {
    compute_load(analyte_obj, analyte, flow_data)
  })

  dplyr::bind_rows(loads)
}

# ============================ RUN LOAD ANALYSIS FUNCTION ============================

warn_if_multiple_sampling_events <- function(wq_df, enable = TRUE) {
  if (!enable || nrow(wq_df) == 0 || !"collected" %in% names(wq_df)) return(invisible(NULL))
  dates <- as.Date(wq_df$collected)
  dates <- dates[is.finite(as.numeric(dates))]
  uniq <- sort(unique(dates))
  if (length(uniq) > 1) {
    counts <- sort(table(dates))
    warning(
      paste0(
        "Multiple WQ sampling events detected in the selected date range. ",
        "The load calculation will average concentrations across these dates. ",
        "Dates (n): ",
        paste(paste0(names(counts), " (", as.integer(counts), ")"), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

run_load_analysis <- function(wq_file, flow_file, start_date, end_date,
                             treatment_filter = "Outflow", user_interval = 4) {

  wq_raw   <- import_wq_data(wq_file)
  print(paste0("WQ Data Imported: ", nrow(wq_raw), " rows"))

  flow_raw <- import_flow_data_flex(flow_file, tz = "America/Denver")
  print(paste0("Flow Data Imported: ", nrow(flow_raw), " rows"))

  wq   <- process_wq_data(wq_raw, treatment_filter)
  print("WQ Data Processed")

  flow <- process_flow_data_flex(flow_raw, tz = "America/Denver")
  print("Flow Data Processed")

  # Accepts either "YYYY-MM-DD" or "YYYY-MM-DD HH:MM" (seconds ok too)
  start_date <- as.POSIXct(trimws(start_date), tz = "America/Denver")
  print(paste0("Start Date: ", start_date))
  end_date <- as.POSIXct(trimws(end_date), tz = "America/Denver")
  print(paste0("End Date: ", end_date))

  wq <- filter(wq, collected >= start_date & collected <= end_date)
  print(paste0("WQ Data Filtered: ", nrow(wq), " rows"))
  if (nrow(wq) == 0) warning("No WQ data found in the selected time frame. Check your treatment filter and date range.", call. = FALSE)
  warn_if_multiple_sampling_events(wq, enable = TRUE)

  flow <- filter(flow, datetime >= start_date & datetime <= end_date)
  print(paste0("Flow Data Filtered: ", nrow(flow), " rows"))
  if (nrow(flow) == 0) warning("No flow data found in the selected time frame. Check your date range.", call. = FALSE)

  flow_aggregated <- aggregate_flow_data(flow, user_interval = user_interval)
  print(paste0("Flow Data Aggregated: ", nrow(flow_aggregated), " rows"))

  flow_sum <- sum(flow_aggregated$volume_L, na.rm = TRUE)
  print(paste0("Total Volume, L: ", flow_sum, " | Gallons: ", flow_sum * 0.264172, " | Acre-ft: ", flow_sum * 8.10714e-7))

  analyte_summary <- compute_analyte_summary(wq)
  analyte_summary <- add_total_nitrogen_to_summary(analyte_summary)

  all_loads <- compute_all_loads(analyte_summary, flow_aggregated)

  print("Excluded non-mass analytes from load calc: EC25, pH")
  if (nrow(all_loads) > 0) {
    print(paste0("Final analytes in load table: ", paste(all_loads$analyte, collapse = ", ")))
  } else {
    print("Final analytes in load table: <none>")
  }

  list(volume = flow_sum, loads = all_loads)
}

# ============================ SUM OBJECTS FXN ============================

sum_load_objects <- function(...) {
  load_objects <- list(...)
  total_volume <- sum(sapply(load_objects, function(obj) obj$volume), na.rm = TRUE)

  original_order <- load_objects[[1]]$loads$analyte

  combined_loads <- bind_rows(lapply(load_objects, function(obj) obj$loads))

  total_loads <- combined_loads %>%
    group_by(analyte) %>%
    summarise(
      load_kg       = sum(load_kg, na.rm = TRUE),
      load_upper_kg = sum(load_upper_kg, na.rm = TRUE),
      load_lower_kg = sum(load_lower_kg, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(analyte, levels = original_order))

  list(volume = total_volume, loads = total_loads)
}

# ================== CONVERT LOAD TO SPATIAL UNITS ============================

# Returns BOTH lb/acre and kg/ha (with bounds), plus acre-ft volume.
convert_load_to_spatial <- function(load_object, acreage) {

  # Conversions
  L_to_acre_ft <- 8.10714e-7
  kg_to_lbs    <- 2.20462
  acre_to_ha   <- 0.40468564224  # 1 acre = 0.40468564224 ha

  area_ha <- acreage * acre_to_ha

  volume_acre_ft <- load_object$volume * L_to_acre_ft

  converted_loads <- load_object$loads %>%
    mutate(
      load_lbs       = load_kg       * kg_to_lbs,
      load_upper_lbs = load_upper_kg * kg_to_lbs,
      load_lower_lbs = load_lower_kg * kg_to_lbs,

      lbs_acre       = load_lbs       / acreage,
      upper_lbs_acre = load_upper_lbs / acreage,
      lower_lbs_acre = load_lower_lbs / acreage,

      kg_ha          = load_kg       / area_ha,
      upper_kg_ha    = load_upper_kg / area_ha,
      lower_kg_ha    = load_lower_kg / area_ha
    ) %>%
    select(
      analyte,
      load_kg, load_upper_kg, load_lower_kg,
      lbs_acre, upper_lbs_acre, lower_lbs_acre,
      kg_ha, upper_kg_ha, lower_kg_ha
    )

  list(
    volume_acre_ft = volume_acre_ft,
    loads = converted_loads
  )
}

# ============================ PLOTLY FLOW TIMESERIES ============================

plot_avg_flow_timeseries_plotly <- function(flow_file,
                                            start_date = NULL,
                                            end_date   = NULL,
                                            user_interval = NULL,   # NULL = no aggregation, else hours
                                            tz = "America/Denver",
                                            flow_units = c("gpm", "cfs", "L/s"),
                                            show_flow_sample_events = TRUE,
                                            sample_event_threshold = 0.5,
                                            wq_file = NULL,
                                            wq_treatment_filter = "Outflow",
                                            show_wq_samples = TRUE,
                                            event_match_window_mins = 60) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install it with install.packages('plotly').", call. = FALSE)
  }
  
  flow_units <- match.arg(flow_units)
  
  # ---- unit conversion (from gpm) ----
  # 1 cfs = 448.831 gpm
  # 1 gpm = 3.78541 L/min = 0.0630901667 L/s
  conv_factor <- switch(
    flow_units,
    "gpm" = 1,
    "cfs" = 1 / 448.831,
    "L/s" = 3.78541 / 60
  )
  y_title <- switch(flow_units, "gpm" = "Flow (gpm)", "cfs" = "Flow (cfs)", "L/s" = "Flow (L/s)")
  y_fmt   <- switch(flow_units, "gpm" = ".2f", "cfs" = ".3f", "L/s" = ".2f")
  
  # Guard against accidentally passing a WQ file as flow_file
  hdr <- tryCatch(read.csv(flow_file, nrows = 1, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(hdr)) {
    nm <- tolower(names(hdr))
    if (all(c("analyte", "result") %in% nm) || ("treatment.name" %in% nm) || ("event.type" %in% nm)) {
      stop(
        paste0(
          "That file looks like a WATER QUALITY table, not an ISCO flow export: ", flow_file, "\n",
          "Pass the flow CSV you use as 'flow_file' in run_load_analysis()."
        ),
        call. = FALSE
      )
    }
  }
  
  if (!is.null(start_date)) start_dt <- as.POSIXct(trimws(start_date), tz = tz) else start_dt <- NULL
  if (!is.null(end_date))   end_dt   <- as.POSIXct(trimws(end_date), tz = tz)   else end_dt   <- NULL
  
  flow_raw <- import_flow_data_flex(flow_file, tz = tz)
  flow     <- process_flow_data_flex(flow_raw, tz = tz)
  
  if (!is.null(start_dt)) flow <- dplyr::filter(flow, datetime >= start_dt)
  if (!is.null(end_dt))   flow <- dplyr::filter(flow, datetime <= end_dt)
  
  if (nrow(flow) == 0) {
    warning("No flow rows found after filtering. Nothing to plot.", call. = FALSE)
    return(invisible(NULL))
  }
  
  # Line layer: optional aggregation
  line_df <- flow
  x_col <- "datetime"
  if (!is.null(user_interval) && is.finite(user_interval) && user_interval > 0) {
    line_df <- aggregate_flow_data(flow, user_interval = user_interval)
    x_col <- "time_block"
  }
  
  # Add plotting y column in desired units (do not overwrite avg_flow_gpm)
  line_df <- line_df %>%
    dplyr::mutate(plot_flow = avg_flow_gpm * conv_factor)
  
  # Build base plotly (flow line)
  p <- plotly::plot_ly(
    data = line_df,
    x = ~ .data[[x_col]],
    y = ~ plot_flow,
    type = "scatter",
    mode = "lines",
    name = paste0("avg_flow_", flow_units),
    hovertemplate = paste0(
      "Time: %{x}",
      "<br>Avg flow (", flow_units, "): %{y:", y_fmt, "}",
      "<extra></extra>"
    )
  )
  
  # ---- Flow sample_event markers (raw flow only) ----
  if (isTRUE(show_flow_sample_events) && ("sample_event" %in% names(flow))) {
    
    ev <- dplyr::filter(flow, is.finite(sample_event) & sample_event >= sample_event_threshold)
    
    if (nrow(ev) > 0) {
      
      flow_valid <- dplyr::filter(flow, is.finite(avg_flow_gpm))
      ev$avg_flow_gpm_marker <- ev$avg_flow_gpm
      
      if (nrow(flow_valid) > 0) {
        
        idx <- vapply(ev$datetime, function(t) {
          j <- which.min(abs(as.numeric(difftime(flow_valid$datetime, t, units = "mins"))))
          if (length(j) == 0) return(NA_integer_)
          if (abs(as.numeric(difftime(flow_valid$datetime[j], t, units = "mins"))) > event_match_window_mins) return(NA_integer_)
          j
        }, FUN.VALUE = integer(1))
        
        ok <- !is.na(idx)
        needs_y <- !is.finite(ev$avg_flow_gpm_marker)
        fill_ok <- ok & needs_y
        ev$avg_flow_gpm_marker[fill_ok] <- flow_valid$avg_flow_gpm[idx[fill_ok]]
      }
      
      ev_plot <- dplyr::filter(ev, is.finite(avg_flow_gpm_marker)) %>%
        dplyr::mutate(
          bottle = as.factor(sample_event),
          plot_flow_marker = avg_flow_gpm_marker * conv_factor
        )
      
      if (nrow(ev_plot) > 0) {
        for (b in sort(unique(ev_plot$bottle))) {
          dfb <- dplyr::filter(ev_plot, bottle == b)
          
          p <- p %>% plotly::add_trace(
            data = dfb,
            x = ~ datetime,
            y = ~ plot_flow_marker,
            type = "scatter",
            mode = "markers",
            name = paste("Bottle", as.character(b)),
            hovertemplate = paste0(
              "Time: %{x}",
              "<br>Avg flow (", flow_units, "): %{y:", y_fmt, "}",
              "<br>Bottle: %{text}",
              "<extra></extra>"
            ),
            text = ~ as.character(sample_event),
            inherit = FALSE
          )
        }
      }
    }
  }
  
  # ---- Optional: overlay WQ sample times as markers ----
  if (!is.null(wq_file) && isTRUE(show_wq_samples)) {
    
    wq_raw <- import_wq_data(wq_file)
    wq <- process_wq_data(wq_raw, wq_treatment_filter)
    
    if (!is.null(start_dt)) wq <- dplyr::filter(wq, collected >= start_dt)
    if (!is.null(end_dt))   wq <- dplyr::filter(wq, collected <= end_dt)
    
    if (nrow(wq) > 0) {
      wq_times <- wq %>%
        dplyr::filter(!is.na(collected)) %>%
        dplyr::distinct(collected) %>%
        dplyr::arrange(collected)
      
      p <- p %>% plotly::add_trace(
        data = wq_times,
        x = ~ collected,
        y = 0,
        type = "scatter",
        mode = "markers",
        name = "WQ collected",
        marker = list(symbol = "x", size = 9),
        hovertemplate = paste(
          "WQ collected: %{x}",
          "<extra></extra>"
        ),
        inherit = FALSE
      )
    }
  }
  
  p <- p %>% plotly::layout(
    title = "Average flow time series",
    xaxis = list(
      title = "Time",
      rangeslider = list(visible = TRUE)
    ),
    yaxis = list(title = y_title),
    hovermode = "x unified"
  )
  
  print(p)
  
  attr(flow, "line_df") <- line_df
  attr(flow, "flow_units") <- flow_units
  invisible(flow)
}
