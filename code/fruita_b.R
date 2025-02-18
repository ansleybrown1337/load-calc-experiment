# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2024 Fruita B Data

B24_event1 <- run_load_analysis(
  wq_file = "./private-data/Fruita B/2024_fruitaB_wq.csv",
  flow_file = "./private-data/Fruita B/2024_fruitaB_flow.csv",
  start_date = "2024-04-13",
  end_date = "2024-04-22",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

B24_event1$volume
B24_event1$loads

B24_event3 <- run_load_analysis(
  wq_file = "./private-data/Fruita B/2024_fruitaB_wq.csv",
  flow_file = "./private-data/Fruita B/2024_fruitaB_flow.csv",
  start_date = "2024-05-21",
  end_date = "2024-5-26",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

B24_event3$volume
B24_event3$loads

total24_B <- sum_load_objects(B24_event1, B24_event3)
total24_B$volume
total24_B$loads

# convert to spatial units
# Assume we have a load object and field size of 50 acres
spatial_loads <- convert_load_to_spatial(total24_B, acreage = 27)

# View total volume in acre-feet
spatial_loads$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads$loads)
