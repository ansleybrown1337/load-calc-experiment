# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2024 Fruita alf Data

alf24_event3 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2024_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2024_fruitaA_flow.csv",
  start_date = "2024-04-13",
  end_date = "2024-04-22",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf24_event3$volume
alf24_event3$loads

alf24_event4 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2024_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2024_fruitaA_flow.csv",
  start_date = "2024-04-13",
  end_date = "2024-04-22",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf24_event4$volume
alf24_event4$loads

alf24_event6 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2024_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2024_fruitaA_flow.csv",
  start_date = "2024-04-13",
  end_date = "2024-04-22",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf24_event6$volume
alf24_event6$loads

alf24_event8 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2024_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2024_fruitaA_flow.csv",
  start_date = "2024-04-13",
  end_date = "2024-04-22",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf24_event8$volume
alf24_event8$loads


total24_alf <- sum_load_objects(alf24_event1, alf24_event3)
total24_alf$volume
total24_alf$loads

# convert to spatial units
# Assume we have a load object and field size of 50 acres
spatial_loads <- convert_load_to_spatial(total24_alf, acreage = 4)

# View total volume in acre-feet
spatial_loads$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads$loads)
