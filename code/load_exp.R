# load-calc-experiment
# Script to calculate water analyte loads and showcase limitations of methods used

# Load required libraries
library(dplyr)
library(tidyr)

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
      collected = as.POSIXct(collected, format="%m/%d/%Y %H:%M", tz="UTC"),
      received = as.POSIXct(received, format="%m/%d/%Y %H:%M", tz="UTC"),
      treatment.name = as.character(treatment.name),
      analyte_abbr = vapply(analyte, function(a) {
        match <- names(analyteAbbr.dict)[sapply(analyteAbbr.dict, function(x) a %in% x)]
        if (length(match) == 0) NA else match
      }, FUN.VALUE = character(1))  # Ensures correct vector length
    ) %>%
    filter(treatment.name == treatment_filter) %>%  # Dynamic filtering
    filter(!is.na(analyte_abbr))  # Remove unmapped analytes
}

## Process Flow Data
process_flow_data <- function(flow_data) {
  # Check how many columns exist
  num_columns <- ncol(flow_data)
  
  # Define column names based on dataset format
  if (num_columns == 8) {
    column_names <- c("datetime", "min_flow_gpm", "min_time", "max_flow_gpm", 
                      "max_time", "avg_flow_gpm", "volume_gal", "sample_event")
  } else if (num_columns == 6) {
    column_names <- c("datetime", "min_flow_gpm", "max_flow_gpm", 
                      "avg_flow_gpm", "volume_gal", "sample_event")
  } else {
    stop("Unexpected number of columns in flow data: ", num_columns)
  }
  
  flow_data <- flow_data %>%
    setNames(column_names) %>%
    mutate(
      datetime = as.POSIXct(datetime, format="%m/%d/%Y %H:%M", tz="UTC"),
      volume_L = if ("volume_gal" %in% names(.)) volume_gal * 3.78541 else NA_real_
    ) %>%
    # Only include time-based columns if they exist
    { if ("min_time" %in% names(.)) 
      mutate(., min_time = as.POSIXct(paste(as.Date(datetime), min_time), 
                                      format="%Y-%m-%d %I:%M:%S %p", tz="UTC")) else . } %>%
    { if ("max_time" %in% names(.)) 
      mutate(., max_time = as.POSIXct(paste(as.Date(datetime), max_time), 
                                      format="%Y-%m-%d %I:%M:%S %p", tz="UTC")) else . }
  
  return(flow_data)
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

## Compute Loads for All Analytes
compute_all_loads <- function(analyte_obj, flow_data) {
  analyte_names <- names(analyte_obj)
  loads <- lapply(analyte_names, function(analyte) compute_load(analyte_obj, analyte, flow_data))
  bind_rows(loads)
}

# ============================ RUN LOAD ANALYSIS FUNCTION ============================

run_load_analysis <- function(wq_file, flow_file, start_date, end_date, treatment_filter = "Outflow", user_interval = 4) {
  
  # ============================ IMPORT & PROCESS DATA ============================
  
  ## Import Data
  wq_raw <- import_wq_data(wq_file)
  flow_raw <- import_flow_data(flow_file)
  
  ## Process Data (now allows user-specified treatment)
  wq <- process_wq_data(wq_raw, treatment_filter)
  flow <- process_flow_data(flow_raw)
  
  ## Filter Data to Date Range
  start_date <- as.POSIXct(start_date, tz = "UTC")
  end_date <- as.POSIXct(end_date, tz = "UTC")
  
  wq <- filter(wq, collected >= start_date & collected <= end_date)
  flow <- filter(flow, datetime >= start_date & datetime <= end_date)
  
  ## Aggregate Flow Data
  flow_aggregated <- aggregate_flow_data(flow, user_interval = user_interval)
  flow_sum <- sum(flow_aggregated$volume_L)
  
  # ============================ ANALYTE PROCESSING ============================
  
  ## Compute Analyte Summary (mean & SD for each analyte)
  analyte_summary <- compute_analyte_summary(wq)
  
  ## Compute Loads for All Analytes
  all_loads <- compute_all_loads(analyte_summary, flow_aggregated)
  
  return(list(volume = flow_sum, loads = all_loads))
}

# ============================ RUN FUNCTION ============================

load_results <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date = "2022-07-08",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_results$volume
load_results$loads
