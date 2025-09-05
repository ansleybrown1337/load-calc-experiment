# load-calc-experiment
# Script to calculate water analyte loads and showcase limitations of methods used

# Load required libraries
library(dplyr)
library(tidyr)
library(lubridate)

# TODO
  # integrate LCS flow data types for import
  # integrate user-based time aggregation in final function

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

# ============================ IMPORT DATA ============================

## Water Quality Data
import_wq_data <- function(file_path) {
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
      # Add support for different formats including '30 May 2024 09:35'
      collected = as.POSIXct(parse_date_time(collected, orders = c("mdY HM", "mdY", "Ymd HMS", "Ymd", "d b Y HM")), tz = "UTC"),
      
      received = as.POSIXct(parse_date_time(received, orders = c("mdY HM", "mdY", "Ymd HMS", "Ymd", "d b Y HM")), tz = "UTC"),
      
      treatment.name = as.character(treatment.name),
      
      # Apply analyte abbreviation mapping
      analyte_abbr = vapply(analyte, function(a) {
        match <- names(analyteAbbr.dict)[sapply(analyteAbbr.dict, function(x) a %in% x)]
        if (length(match) == 0) NA else match
      }, FUN.VALUE = character(1))
    ) %>%
    # Adjust treatment filter to check both 'treatment.name' and 'event.type'
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
  
  flow_data <- flow_data %>%
    setNames(column_names) %>%
    mutate(
      # parse timestamp
      datetime = as.POSIXct(datetime, format = "%m/%d/%Y %H:%M", tz = "UTC"),
      # clean and coerce numerics
      across(
        intersect(c("min_flow_gpm","max_flow_gpm","avg_flow_gpm","volume_gal","sample_event"), names(.)),
        ~ suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(.)))))
      ),
      # now safe to compute liters
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
  
  flow_data
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
    # Replace Inf values with NA to prevent issues
    mutate(
      min_flow_gpm = ifelse(is.infinite(min_flow_gpm), NA_real_, min_flow_gpm),
      max_flow_gpm = ifelse(is.infinite(max_flow_gpm), NA_real_, max_flow_gpm)
    ) %>%
    # Remove NA columns only if they weren’t originally present
    select(where(~ !all(is.na(.))))
}



# ============================ COMPUTE ANALYTE SUMMARY ============================

compute_analyte_summary <- function(wq_data) {
  wq_data %>%
    group_by(analyte_abbr) %>%
    summarise(
      avg_result = mean(result, na.rm = TRUE),
      sd_result = sd(result, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    split(.$analyte_abbr) %>%
    lapply(function(df) list(avg = df$avg_result, sd = df$sd_result))
}

# ============================ LOAD CALCULATION ============================

## Compute Load for a Given Analyte
compute_load <- function(analyte_obj, analyte_name, flow_data) {
  if (!analyte_name %in% names(analyte_obj)) {
    return(NA)  # Return NA if the analyte is not found
  }
  
  avg_conc <- analyte_obj[[analyte_name]]$avg
  sd_conc <- analyte_obj[[analyte_name]]$sd
  total_volume <- sum(flow_data$volume_L)
  
  # Load calculations
  load_kg <- avg_conc * total_volume / 1e6  # Convert mg to kg
  load_kg_upper <- (avg_conc + sd_conc) * total_volume / 1e6
  load_kg_lower <- (avg_conc - sd_conc) * total_volume / 1e6
  
  return(data.frame(
    analyte = analyte_name,
    load_kg = load_kg,
    load_upper_kg = load_kg_upper,
    load_lower_kg = load_kg_lower
  ))
}

## Compute Loads for All Analytes (with exclusions for ec&ph + status message)
compute_all_loads <- function(analyte_obj, flow_data, drop_analytes = c("EC25", "pH")) {
  # Keep only analytes that should have mass-based loads
  analyte_names <- setdiff(names(analyte_obj), drop_analytes)
  
  # Status message if any analytes are dropped
  if (length(drop_analytes) > 0) {
    message("Note: The following analytes were excluded from load calculations because they are not mass-based: ",
            paste(drop_analytes, collapse = ", "))
  }
  
  # If nothing left to compute, return an empty-but-correctly-shaped frame
  if (length(analyte_names) == 0) {
    return(dplyr::tibble(
      analyte = character(),
      load_kg = double(),
      load_upper_kg = double(),
      load_lower_kg = double()
    ))
  }
  
  # Compute loads
  loads <- lapply(analyte_names, function(analyte) {
    compute_load(analyte_obj, analyte, flow_data)
  })
  
  dplyr::bind_rows(loads)
}



# ============================ RUN LOAD ANALYSIS FUNCTION ============================

run_load_analysis <- function(wq_file, flow_file, start_date, end_date, treatment_filter = "Outflow", user_interval = 4) {
  
  # ============================ IMPORT & PROCESS DATA ============================
  
  ## Import Data
  wq_raw <- import_wq_data(wq_file)
  print(paste0("WQ Data Imported: ", nrow(wq_raw), " rows"))
  flow_raw <- import_flow_data(flow_file)
  print(paste0("Flow Data Imported: ", nrow(flow_raw), " rows"))
  
  ## Process Data (now allows user-specified treatment)
  wq <- process_wq_data(wq_raw, treatment_filter)
  print("WQ Data Processed")
  flow <- process_flow_data(flow_raw)
  print("Flow Data Processed")
  
  ## Filter Data to Date Range
  ### aj note 10 feb 2024: doesn't filter well if NAs present in dates
  start_date <- as.POSIXct(start_date, tz = "UTC")
  print(paste0("Start Date: ", start_date))
  end_date <- as.POSIXct(end_date, tz = "UTC")
  print(paste0("End Date: ", end_date))
  
  wq <- filter(wq, collected >= start_date & collected <= end_date)
  print(paste0("WQ Data Filtered: ", nrow(wq), " rows"))
  if (nrow(wq) == 0) {
    warning("No WQ data found in the selected time frame. Check your treatment filter and date range.")
  }
  flow <- filter(flow, datetime >= start_date & datetime <= end_date)
  print(paste0("Flow Data Filtered: ", nrow(flow), " rows"))
  if (nrow(flow) == 0) {
    warning("No flow data found in the selected time frame. Check your treatment filter and date range.")
  }
  
  ## Aggregate Flow Data
  flow_aggregated <- aggregate_flow_data(flow, user_interval = user_interval)
  print(paste0("Flow Data Aggregated: ", nrow(flow_aggregated), " rows"))
  flow_sum <- sum(flow_aggregated$volume_L)
  print(paste0("Total Volume, L: ", flow_sum, " | Gallons: ", flow_sum * 0.264172))
  
  # ============================ ANALYTE PROCESSING ============================
  
  ## Compute Analyte Summary (mean & SD for each analyte)
  analyte_summary <- compute_analyte_summary(wq)
  
  ## Compute Loads for All Analytes
  all_loads <- compute_all_loads(analyte_summary, flow_aggregated)
  
  return(list(volume = flow_sum, loads = all_loads))
}

# ============================ RUN FUNCTION EXAMPLE ============================

# load_results <- run_load_analysis(
#   wq_file = "./data/2022uym_wq.csv",
#   flow_file = "./data/2022uym_flow.csv",
#   start_date = "2022-06-12",
#   end_date = "2022-07-08",
#   treatment_filter = "Outflow",  # User can now specify treatment!
#   user_interval = 4
# )
# 
# # Access results
# load_results$volume
# load_results$loads

# ============================ SUM OBJECTS FXN ============================

sum_load_objects <- function(...) {
  # Capture all load objects passed to the function
  load_objects <- list(...)
  
  # Sum total volume
  total_volume <- sum(sapply(load_objects, function(obj) obj$volume), na.rm = TRUE)
  
  # Extract original analyte order from the first load object
  original_order <- load_objects[[1]]$loads$analyte
  
  # Combine all loads into a single data frame
  combined_loads <- bind_rows(lapply(load_objects, function(obj) obj$loads))
  
  # Sum loads by analyte and preserve original order
  total_loads <- combined_loads %>%
    group_by(analyte) %>%
    summarise(
      load_kg = sum(load_kg, na.rm = TRUE),
      load_upper_kg = sum(load_upper_kg, na.rm = TRUE),
      load_lower_kg = sum(load_lower_kg, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Preserve the original order of analytes
    arrange(factor(analyte, levels = original_order))
  
  # Return the summed object with the same structure
  summed_object <- list(
    volume = total_volume,
    loads = total_loads
  )
  
  return(summed_object)
}

# ================== CONVERT LOAD TO SPATIAL UNITS ============================

convert_load_to_spatial <- function(load_object, acreage) {
  
  # Conversion factors
  
  L_to_acre_ft <- 8.10714e-7  # 1 liter = 8.10714e-7 acre-feet
  kg_to_lbs <- 2.20462      # 1 kg = 2.20462 lbs
  g_to_lbs <- 0.00220462    # 1 g = 0.00220462 lbs
  
  # Convert volume to acre-feet
  
  volume_acre_ft <- load_object$volume * L_to_acre_ft
  
  # Convert loads to lbs, kg, and lbs/acre
  
  converted_loads <- load_object$loads %>%
    mutate(
      load_lbs = ifelse(
        analyte %in% c("Fe (g not kg)", "Se (g not kg)"),  # If analyte is in grams
        load_kg * g_to_lbs,                                  # Convert g → lbs
        load_kg * kg_to_lbs                                  # Convert kg → lbs
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
      lbs_acre = load_lbs / acreage,
      upper_lbs_acre = load_upper_lbs / acreage,
      lower_lbs_acre = load_lower_lbs / acreage
    ) %>%
    select(analyte,load_lbs, load_upper_lbs, load_lower_lbs,
           lbs_acre, upper_lbs_acre, lower_lbs_acre)
  
  # Return converted object with the same structure
  
  return(list(
    volume_acre_ft = volume_acre_ft,
    loads = converted_loads
  ))
}

