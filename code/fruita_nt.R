# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2024 Fruita NT Data

# event 1 - grab only, no loads

nt24_event2 <- run_load_analysis(
  wq_file = "./private-data/Fruita NT/2024_fruitaNT_wq.csv",
  flow_file = "./private-data/Fruita NT/2024_fruitaNT_flow.csv",
  start_date = "2024-04-15",
  end_date = "2024-07-01",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

nt24_event2$volume
nt24_event2$loads

nt24_event5 <- run_load_analysis(
  wq_file = "./private-data/Fruita NT/2024_fruitaNT_wq.csv",
  flow_file = "./private-data/Fruita NT/2024_fruitaNT_flow.csv",
  start_date = "2024-09-15",
  end_date = "2024-10-13",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

nt24_event5$volume
nt24_event5$loads

total24_nt <- sum_load_objects(nt24_event2, nt24_event5)
total24_nt$volume
total24_nt$loads

# convert to spatial units
# Assume we have a load object and field size of 50 acres
spatial_loads <- convert_load_to_spatial(total24_nt, acreage = 20)

# View total volume in acre-feet
spatial_loads$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads$loads)
