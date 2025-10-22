# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/Jake/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2025 Fruita alf Data

alf25_event1 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-04-14",
  end_date = "2025-04-16",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event1$volume
alf25_event1$loads

alf25_event2 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-05-10",
  end_date = "2025-05-14",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event2$volume
alf25_event2$loads

alf25_event3 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-06-15",
  end_date = "2025-06-18",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event3$volume
alf25_event3$loads

alf25_event4 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-06-25",
  end_date = "2025-06-28",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event4$volume
alf25_event4$loads

alf25_event5 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-07-12",
  end_date = "2025-07-16",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event5$volume
alf25_event5$loads

alf25_event6 <- run_load_analysis(
  wq_file = "./private-data/Fruita A/2025_fruitaA_wq.csv",
  flow_file = "./private-data/Fruita A/2025_fruitaA_flow.csv",
  start_date = "2025-08-26",
  end_date = "2025-08-30",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

alf25_event6$volume
alf25_event6$loads


total25_alf <- sum_load_objects(alf25_event1, alf25_event2, alf25_event3, alf25_event4, alf25_event5, alf25_event6)
total25_alf$volume
total25_alf$loads

# convert to spatial units
# Assume we have a load object and field size of 50 acres
spatial_loads <- convert_load_to_spatial(total25_alf, acreage = 4)

# View total volume in acre-feet
spatial_loads$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads$loads)
write.csv(spatial_loads, "./code/spatial_loads.csv")
