# load-calc-experiment
# Script to calculate water analyte loads and showcase limitations of methods used

# Load required libraries
library(dplyr)
library(tidyr)
library(lubridate)

# TODO
#   - integrate LCS flow data types for import
#   - integrate user-based time aggregation in final function

# ============================ ANALYTE ABBREVIATION DICTIONARY ============================

analyteAbbr.dict <- list(
  "TKN"   = c("Nitrogen, Total Kjeldahl"),
  "NO2_N" = c("Nitrogen, Nitrite  (As N)", "NITRITE AS N"),
  "PO4_P" = c("Phosphorus, Total Orthophosphate (as P)", "ORTHOPHOSPHATE AS P"),
  "TP"    = c("Phosphorus, Total (As P)", "TOTAL PHOSPHORUS"),
  "TDS"   = c("Total Dissolved Solids (Residue, Filterable)"),
  "NO3_N" = c("Nitrogen, Nitrate (As N)", "NITRATE AS N"),
  "TSS"   = c("Suspended Solids (Residue, Non-Filterable)", "TSS"),
  "Fe (g not kg)"    = c("Iron, Total"),
  "Se (g not kg)"    = c("Selenium, Total", "SELENIUM"),
  "pH"    = c("pH", "POTENTIAL HYDROGEN"),
  "EC25"  = c("Specific Conductance", "ELECTRICAL CONDUCTIVITY")
)

# ============================ SMALL HELPERS ============================

# Clean numeric-ish columns (removes commas/whitespace and coerces to numeric)
clean_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))
}

# Parse censored values / strings in WQ result column
# NOTE: Per user request, NDs and values like "<0.01" are treated as ZERO.
parse_result_to_numeric_nd_as_zero <- function(x) {
  x_chr <- trimws(as.character(x))
  is_lt <- grepl("^<", x_chr)
  # strip leading "<" then coerce
  mdls  <- suppressWarnings(as.numeric(sub("^<", "", x_chr)))
  nums  <- suppressWarnings(as.numeric(x_chr))
  out   <- ifelse(is_lt, 0, nums)      # <MDL -> 0
  # common ND tokens -> 0
  out[grepl("^\\s*(ND|nd|non[- ]?detect(ed)?|n/?d)\\s*$", x_chr)] <- 0
  out[!is.finite(out)] <- NA_real_
  out
}

# ============================ IMPORT DATA ============================

## Water Quality Data
import_wq_data <- function(file_path) {
  # keep strings as factors to avoid breaking current expectations
  read.csv(file_path, stringsAsFactors = TRUE)
}

## Flow Data
import_flow_data <- function(file_path) {
  read.csv(file_path, skip = 7, header = FALSE, stringsAsFactors = FALSE)
}

# ============================ PROCESS DATA ============================

## Process Water Quality Data (with treatment filter)
process_wq_data <- function(wq_data, treatment_filter) {
  wq_data %>%
    mutate(
      # Robust timestamp parsing, e.g. '30 May 2024 09:35'
      collected = as.POSIXct(parse_date_time(
        collected, orders = c("mdY HM","mdY","Ymd HMS","Ymd","d b Y HM")
      ), tz = "UTC"),
      received = as.POSIXct(parse_date_time(
        received, orders = c("mdY HM","mdY","Ymd HMS","Ymd","d b Y HM")
      ), tz = "UTC"),
      
      treatment.name = as.character(treatment.name),
      
      # Map to analyte abbreviations
      analyte_abbr = vapply(analyte, function(a) {
        match <- names(analyteAbbr.dict)[sapply(analyteAbbr.dict, function(x) a %in% x)]
        if (length(match) == 0) NA_character_ else match
      }, FUN.VALUE = character(1)),
      
      # Convert results to numeric with ND handling
      # NOTE: ND and "<MDL" values are treated as ZERO by request
      result = parse_result_to_numeric_nd_as_zero(result)
    ) %>%
    {
      if (all(is.na(.$treatment.name))) {
        stop("Warning: All values in 'treatment.name' are NA. Please verify your dataset and ensure 'treatment.name' is correctly populated.")
      } else {
        filter(., treatment.name == treatment_filter | event.type == treatment_filter)
      }
    } %>%
    {
      if (any(is.na(.$analyte_abbr))) {
        stop("Warning: Some analytes could not be mapped. Please check for unmapped analytes and update 'analyteAbbr.dict' accordingly.")
      } else {
        .
      }
    }
}

## Process Flow Data
process_flow_data <- function(flow_data) {
  num_columns <- ncol(flow_data)
  
  if (num_columns == 8) {
    column_names <- c("datetime","min_flow_gpm","min_time","max_flow_gpm",
                      "max_time","avg_flow_gpm","volume_gal","sample_event")
  } else if (num_columns == 6) {
    column_names <- c("datetime","min_flow_gpm","max_flow_gpm",
                      "avg_flow_gpm","volume_gal","sample_event")
  } else {
    stop("Unexpected number of columns in flow data: ", num_columns)
  }
  
  flow_data %>%
    setNames(column_names) %>%
    mutate(
      # parse timestamp
      datetime = as.POSIXct(datetime, format = "%m/%d/%Y %H:%M", tz = "UTC"),
      # clean and coerce numerics
      across(
        intersect(c("min_flow_gpm","max_flow_gpm","avg_flow_gpm","volume_gal","sample_event"), names(.)),
        clean_numeric
      ),
      # liters from gallons if available
      volume_L = if ("volume_gal" %in% names(.)) volume_gal * 3.78541 else NA_real_
    ) %>%
    { if ("min_time" %in% names(.))
      mutate(., min_time = as.POSIXct(paste(as.Date(datetime), min_time),
                                      format = "%Y-%m-%d %I:%M:%S %p", tz = "UTC"))
      else . } %>%
    { if ("max_time" %in% names(.))
      mutate(., max_time = as.POSIXct(paste(as.Date(datetime), max_time),
                                      format = "%Y-%m-%d %I:%M:%S %p", tz = "UTC"))
      else . }
}

# ============================ AGGREGATE FLOW DATA ============================

## Aggregate Flow Data (Default: 4-hour interval)
aggregate_flow_data <- function(flow_data, user_interval = 4) {
  flow_data %>%
    mutate(
      time_block = as.POSIXct(
        floor(as.numeric(datetime) / (user_interval * 3600)) * (user_interval * 3600),
        origin = "1970-01-01", tz = "UTC"
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
      sample_event = sum(sample_event, na.rm = TRUE),
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

compute_analyte_summary <- function(wq_data) {
  wq_data %>%
    group_by(analyte_abbr) %>%
    summarise(
      avg_result = mean(result, na.rm = TRUE),
      sd_result  = sd(result, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    split(.$analyte_abbr) %>%
    lapply(function(df) list(avg = df$avg_result, sd = df$sd_result))
}

# -------- PROTOTYPE (not used): Flow-Weighted Mean Concentrations (FWMC) --------
# This is a prototype; not wired into run_load_analysis() to avoid any change in behavior.
align_wq_to_flow_prototype <- function(wq, flow_agg, max_diff_mins = 120) {
  # nearest join within tolerance
  idx <- sapply(wq$collected, function(t) {
    if (is.na(t)) return(NA_integer_)
    which.min(abs(difftime(flow_agg$time_block, t, units = "mins")))
  })
  nearest <- ifelse(is.na(idx), NA, flow_agg$time_block[pmin(pmax(idx, 1L), nrow(flow_agg))])
  ok <- !is.na(nearest) &
    abs(difftime(nearest, wq$collected, units = "mins")) <= max_diff_mins
  wq$time_block <- as.POSIXct(NA, tz = "UTC", origin = "1970-01-01")
  wq$time_block[ok] <- nearest[ok]
  dplyr::filter(wq, !is.na(time_block))
}

compute_analyte_summary_fwmc_prototype <- function(wq_aligned, flow_agg) {
  dplyr::left_join(wq_aligned, flow_agg[, c("time_block","volume_L")], by = "time_block") %>%
    group_by(analyte_abbr) %>%
    summarise(
      avg = sum(result * volume_L, na.rm = TRUE) / sum(volume_L, na.rm = TRUE),
      sd  = sd(result, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    split(.$analyte_abbr) %>%
    lapply(function(df) list(avg = df$avg, sd = df$sd))
}
# -------- END PROTOTYPE --------

# ============================ LOAD CALCULATION ============================

## Compute Load for a Given Analyte
compute_load <- function(analyte_obj, analyte_name, flow_data) {
  if (!analyte_name %in% names(analyte_obj)) return(NA)
  
  avg_conc     <- analyte_obj[[analyte_name]]$avg
  sd_conc      <- analyte_obj[[analyte_name]]$sd
  total_volume <- sum(flow_data$volume_L, na.rm = TRUE)
  
  if (!is.finite(total_volume) || total_volume <= 0) {
    warning("Total flow volume is non-finite or <= 0; loads returned as NA for analyte: ", analyte_name)
    return(data.frame(
      analyte = analyte_name,
      load_kg = NA_real_, load_upper_kg = NA_real_, load_lower_kg = NA_real_
    ))
  }
  
  upper <- (avg_conc + sd_conc)
  lower <- pmax(0, avg_conc - sd_conc)  # clamp lower bound to zero
  
  data.frame(
    analyte        = analyte_name,
    load_kg        = avg_conc * total_volume / 1e6,
    load_upper_kg  = upper    * total_volume / 1e6,
    load_lower_kg  = lower    * total_volume / 1e6
  )
}

## Compute Loads for All Analytes (exclude EC25/pH + status message)
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

# Warn if multiple sampling dates are present in the filtered WQ
warn_if_multiple_sampling_events <- function(wq_df, enable = TRUE) {
  if (!enable || nrow(wq_df) == 0 || !"collected" %in% names(wq_df)) return(invisible(NULL))
  # collapse to date (not datetime) so dup samples on same day don't count as separate events
  dates <- as.Date(wq_df$collected)
  dates <- dates[is.finite(as.numeric(dates))]  # drop NA dates quietly
  uniq <- sort(unique(dates))
  if (length(uniq) > 1) {
    # counts by date for context
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


run_load_analysis <- function(wq_file, flow_file, start_date, end_date, treatment_filter = "Outflow", user_interval = 4) {
  
  # ============================ IMPORT & PROCESS DATA ============================
  
  ## Import Data
  wq_raw   <- import_wq_data(wq_file)
  print(paste0("WQ Data Imported: ", nrow(wq_raw), " rows"))
  flow_raw <- import_flow_data(flow_file)
  print(paste0("Flow Data Imported: ", nrow(flow_raw), " rows"))
  
  ## Process Data (now allows user-specified treatment)
  wq   <- process_wq_data(wq_raw, treatment_filter)
  print("WQ Data Processed")
  flow <- process_flow_data(flow_raw)
  print("Flow Data Processed")
  
  ## Filter Data to Date Range
  # aj note 10 feb 2024: doesn't filter well if NAs present in dates
  start_date <- as.POSIXct(start_date, tz = "UTC")
  print(paste0("Start Date: ", start_date))
  end_date <- as.POSIXct(end_date, tz = "UTC")
  print(paste0("End Date: ", end_date))
  
  wq   <- filter(wq, collected >= start_date & collected <= end_date)
  print(paste0("WQ Data Filtered: ", nrow(wq), " rows"))
  if (nrow(wq) == 0) {
    warning("No WQ data found in the selected time frame. Check your treatment filter and date range.")
  }
  # Inform user if multiple sampling dates will be averaged
  warn_if_multiple_sampling_events(wq, enable = TRUE)
  
  
  flow <- filter(flow, datetime >= start_date & datetime <= end_date)
  print(paste0("Flow Data Filtered: ", nrow(flow), " rows"))
  if (nrow(flow) == 0) {
    warning("No flow data found in the selected time frame. Check your treatment filter and date range.")
  }
  
  ## Aggregate Flow Data
  flow_aggregated <- aggregate_flow_data(flow, user_interval = user_interval)
  print(paste0("Flow Data Aggregated: ", nrow(flow_aggregated), " rows"))
  flow_sum <- sum(flow_aggregated$volume_L, na.rm = TRUE)
  print(paste0("Total Volume, L: ", flow_sum, " | Gallons: ", flow_sum * 0.264172))
  
  # ============================ ANALYTE PROCESSING ============================
  
  ## Compute Analyte Summary (mean & SD for each analyte)
  analyte_summary <- compute_analyte_summary(wq)
  
  ## Compute Loads for All Analytes (excludes EC25/pH internally)
  all_loads <- compute_all_loads(analyte_summary, flow_aggregated)
  
  ## Always print a status line about exclusions + final analytes shown
  print("Excluded non-mass analytes from load calc: EC25, pH")
  if (nrow(all_loads) > 0) {
    print(paste0("Final analytes in load table: ", paste(all_loads$analyte, collapse = ", ")))
  } else {
    print("Final analytes in load table: <none>")
  }
  
  return(list(volume = flow_sum, loads = all_loads))
}

# ============================ SUM OBJECTS FXN ============================

sum_load_objects <- function(...) {
  load_objects <- list(...)
  total_volume <- sum(sapply(load_objects, function(obj) obj$volume), na.rm = TRUE)
  
  # Extract original analyte order from the first load object
  original_order <- load_objects[[1]]$loads$analyte
  
  # Combine & sum
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

convert_load_to_spatial <- function(load_object, acreage) {
  
  # Conversion factors
  L_to_acre_ft <- 8.10714e-7  # 1 liter = 8.10714e-7 acre-feet
  kg_to_lbs    <- 2.20462     # 1 kg = 2.20462 lbs
  g_to_lbs     <- 0.00220462  # 1 g = 0.00220462 lbs
  
  # Convert volume to acre-feet
  volume_acre_ft <- load_object$volume * L_to_acre_ft
  
  # Convert loads to lbs and lbs/acre
  converted_loads <- load_object$loads %>%
    mutate(
      load_lbs = ifelse(
        analyte %in% c("Fe (g not kg)", "Se (g not kg)"),
        load_kg * g_to_lbs,  # g -> lbs
        load_kg * kg_to_lbs  # kg -> lbs
      ),
      load_upper_lbs = ifelse(
        analyte %in% c("Fe (g not kg)", "Se (g not kg)"),
        load_upper_kg * g_to_lbs,
        load_upper_kg * kg_to_lbs
      ),
      load_lower_lbs = ifelse(
        analyte %in% c("Fe (g not kg)", "Se (g not kg)"),
        load_lower_kg * g_to_lbs,
        load_lower_kg * kg_to_lbs
      ),
      lbs_acre       = load_lbs       / acreage,
      upper_lbs_acre = load_upper_lbs / acreage,
      lower_lbs_acre = load_lower_lbs / acreage
    ) %>%
    select(analyte, load_lbs, load_upper_lbs, load_lower_lbs,
           lbs_acre, upper_lbs_acre, lower_lbs_acre)
  
  list(
    volume_acre_ft = volume_acre_ft,
    loads = converted_loads
  )
}
