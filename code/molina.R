# working file to use the load_exp.R functions on real data

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2022 SC Data

mol22 <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date = "2022-07-08",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

mol22$volume
mol22$loads

## 2023 SC Data


total23_uym <- sum_load_objects(load_uym23_event1, load_uym23_event2, load_uym23_event3)
total23_uym$volume
total23_uym$loads

## 2024 SC Data

