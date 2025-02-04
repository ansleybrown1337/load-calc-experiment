# load-calc-experiment

# Script to calculate water analyte loads and showcase limitations of methods used

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Import Data
## Water quality data
import_wq_data <- function(file_path) {
  read.csv(file_path, stringsAsFactors = TRUE)
}

## Flow data
import_flow_data <- function(file_path) {
  read.csv(file_path, skip = 7, header = FALSE, stringsAsFactors = FALSE)
}

# Process Data
## Placeholder for water quality processing
process_wq_data <- function(wq_data) {
  wq_data %>%
    mutate(
      collected = as.POSIXct(collected, format="%m/%d/%Y %H:%M", tz="UTC"),
      received = as.POSIXct(received, format="%m/%d/%Y %H:%M", tz="UTC"),
      treatment.name = as.character(treatment.name)  # Convert factor to character
    ) %>%
    filter(treatment.name == "Outflow")  # Now filter safely
}


## Process flow data
process_flow_data <- function(flow_data) {
  column_names <- c("datetime", "min_flow_gpm", "min_time", "max_flow_gpm", 
                    "max_time", "avg_flow_gpm", "volume_gal", "sample_event")
  
  flow_data %>%
    setNames(column_names) %>%
    mutate(
      datetime = as.POSIXct(datetime, format="%m/%d/%Y %H:%M", tz="UTC"),
      date_only = as.Date(datetime),  # Extract just the date part
      min_time = as.POSIXct(paste(date_only, min_time), format="%Y-%m-%d %I:%M:%S %p", tz="UTC"),
      max_time = as.POSIXct(paste(date_only, max_time), format="%Y-%m-%d %I:%M:%S %p", tz="UTC"),
      volume_L = volume_gal * 3.78541  # Convert volume to liters
    ) %>%
    select(-date_only) # Remove intermediate column
}

# Aggregate flow data on sample interval timescale (1 hour default, 4 hours used for this site)
avg_flow_rates <- function(flow_data, user_interval = 1) {
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
      volume_gal = sum(volume_gal, na.rm = TRUE),
      volume_L = sum(volume_gal, na.rm = TRUE),
      sample_event = sum(sample_event, na.rm = TRUE),
      .groups = "drop"
    )
}

# Obtain Average concentration of analytes
compute_analyte_summary <- function(wq_data) {
  wq_data %>%
    group_by(analyte) %>%
    summarise(
      avg_result = mean(result, na.rm = TRUE),
      sd_result = sd(result, na.rm = TRUE),  # Compute standard deviation
      .groups = "drop"
    ) %>%
    split(.$analyte) %>%  # Create a list of lists, indexed by analyte names
    lapply(function(df) list(avg = df$avg_result, sd = df$sd_result))
}

# Load method 1: using analyte mean and total flow volume
load_method_1 <- function(analyte_obj, analyte_index, flow_data) {
  # calculate load
  load_mg <- analyte_obj[[analyte_index]]$avg * sum(flow_data$volume_L)
  # calculate load upper and lower bounds with standard deviation
  load_mg_upper <- (analyte_obj[[analyte_index]]$avg + analyte_obj[[analyte_index]]$sd) * sum(flow_data$volume_L)
  load_mg_lower <- (analyte_obj[[analyte_index]]$avg - analyte_obj[[analyte_index]]$sd) * sum(flow_data$volume_L)
  # convert everything to kg
  load_kg <- load_mg / 1000 / 1000
  load_kg_upper <- load_mg_upper / 1000 / 1000
  load_kg_lower <- load_mg_lower / 1000 / 1000
  # print out analyte chosen
  print(avg_analyte_results[analyte_index])
  # print load with upper and lower bounds
  print(paste0("Load (kg): ", load_kg))
  print(paste0("Load upper bound (kg): ", load_kg_upper))
  print(paste0("Load lower bound (kg): ", load_kg_lower))
  load_kg
}

########### Execute functions ###########
## Import raw data
wq_raw <- import_wq_data("./data/2022uym_wq.csv")
flow_raw <- import_flow_data("./data/2022uym_flow.csv")

## Process data
wq <- process_wq_data(wq_raw)
flow <- process_flow_data(flow_raw)

## Cut data to date range of interest
start_date <- as.POSIXct("2022-06-12", tz = "UTC")
end_date <- as.POSIXct("2022-07-08", tz = "UTC")

wq <- wq %>%
  filter(collected >= start_date & collected <= end_date)

flow <- flow %>%
  filter(datetime >= start_date & datetime <= end_date)

## aggregate flow data
flow_aggregated <- avg_flow_rates(flow, user_interval = 4)

## Compute average analyte results
analyte_summary <- compute_analyte_summary(wq)
print(analyte_summary)

## Load calc 1: total P in kg
load_P_m1 <- load_method_1(analyte_summary, 8, flow_aggregated)

## Load calc 1: total N02 in kg
load_N_m1 <- load_method_1(analyte_summary, 2, flow_aggregated)

## Load calc 1: total TSS in kg
load_TSS_m1 <- load_method_1(analyte_summary, 7, flow_aggregated)
