# load-calc-experiment
# Script to calculate water analyte loads and showcase limitations of methods used

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

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

## Process Water Quality Data
process_wq_data <- function(wq_data) {
  wq_data %>%
    mutate(
      collected = as.POSIXct(collected, format="%m/%d/%Y %H:%M", tz="UTC"),
      received = as.POSIXct(received, format="%m/%d/%Y %H:%M", tz="UTC"),
      treatment.name = as.character(treatment.name)
    ) %>%
    filter(treatment.name == "Outflow")
}

## Process Flow Data
process_flow_data <- function(flow_data) {
  column_names <- c("datetime", "min_flow_gpm", "min_time", "max_flow_gpm", 
                    "max_time", "avg_flow_gpm", "volume_gal", "sample_event")
  
  flow_data %>%
    setNames(column_names) %>%
    mutate(
      datetime = as.POSIXct(datetime, format="%m/%d/%Y %H:%M", tz="UTC"),
      date_only = as.Date(datetime),
      min_time = as.POSIXct(paste(date_only, min_time), format="%Y-%m-%d %I:%M:%S %p", tz="UTC"),
      max_time = as.POSIXct(paste(date_only, max_time), format="%Y-%m-%d %I:%M:%S %p", tz="UTC"),
      volume_L = volume_gal * 3.78541  # Convert gallons to liters
    ) %>%
    select(-date_only)
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
      min_flow_gpm = min(min_flow_gpm, na.rm = TRUE),
      max_flow_gpm = max(max_flow_gpm, na.rm = TRUE),
      avg_flow_gpm = mean(avg_flow_gpm, na.rm = TRUE),
      volume_L = sum(volume_L, na.rm = TRUE),
      sample_event = sum(sample_event, na.rm = TRUE),
      .groups = "drop"
    )
}

# ============================ COMPUTE ANALYTE SUMMARY ============================

compute_analyte_summary <- function(wq_data) {
  wq_data %>%
    group_by(analyte) %>%
    summarise(
      avg_result = mean(result, na.rm = TRUE),
      sd_result = sd(result, na.rm = TRUE),  # Standard deviation for uncertainty
      .groups = "drop"
    ) %>%
    split(.$analyte) %>%
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
  bind_rows(loads)  # Combine into a single dataframe
}

# ============================ PLOTTING FUNCTION ============================

## Create Faceted Bar Graph of Loads for Each Analyte
# mostly placeholder; not working and not including treatment colors yet
plot_analyte_loads <- function(load_data) {
  ggplot(load_data, aes(x = analyte, y = load_kg, fill = analyte)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = load_lower_kg, ymax = load_upper_kg), width = 0.2, color = "black") +
    labs(title = "Analyte Load Estimates", x = "Analyte", y = "Load (kg)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~analyte, scales = "free_y")  # Facet by analyte for clarity
}

# ============================ EXECUTION ============================

## Import Data
wq_raw <- import_wq_data("./data/2022uym_wq.csv")
flow_raw <- import_flow_data("./data/2022uym_flow.csv")

## Process Data
wq <- process_wq_data(wq_raw)
flow <- process_flow_data(flow_raw)

## Filter Data to Date Range
start_date <- as.POSIXct("2022-06-12", tz = "UTC")
end_date <- as.POSIXct("2022-07-08", tz = "UTC")

wq <- filter(wq, collected >= start_date & collected <= end_date)
flow <- filter(flow, datetime >= start_date & datetime <= end_date)

## Aggregate Flow Data
flow_aggregated <- aggregate_flow_data(flow, user_interval = 4)

## Compute Analyte Summary
analyte_summary <- compute_analyte_summary(wq)

## Compute Loads for All Analytes
all_loads <- compute_all_loads(analyte_summary, flow_aggregated)

## Print Results
print(all_loads)
