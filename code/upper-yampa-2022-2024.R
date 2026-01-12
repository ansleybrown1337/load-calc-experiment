# working file to use the load_exp.R functions on real data

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2022 UYM Data

load_uym22 <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date = "2022-07-08",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_uym22$volume
load_uym22$loads

## 2023 UYM Data

load_uym23_event1 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2023uym_wq.csv",
  flow_file = "./private-data/Yampa/2023uym_flow.csv",
  start_date = "2023-06-01",
  end_date = "2023-06-14",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_uym23_event1$volume
load_uym23_event1$loads

load_uym23_event2 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2023uym_wq.csv",
  flow_file = "./private-data/Yampa/2023uym_flow.csv",
  start_date = "2023-06-15",
  end_date = "2023-07-05",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_uym23_event2$volume
load_uym23_event2$loads

load_uym23_event3 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2023uym_wq.csv",
  flow_file = "./private-data/Yampa/2023uym_flow.csv",
  start_date = "2023-07-06",
  end_date = "2023-08-15",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_uym23_event3$volume
load_uym23_event3$loads

total23_uym <- sum_load_objects(load_uym23_event1, load_uym23_event2, load_uym23_event3)
total23_uym$volume
total23_uym$loads

## 2024 UYM Data
load_uym24_event1 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2024uym_wq.csv",
  flow_file = "./private-data/Yampa/2024uym_flow.csv",
  start_date = "2024-05-01",
  end_date = "2024-06-08",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)
load_uym24_event1$volume
load_uym24_event1$loads

load_uym24_event2 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2024uym_wq.csv",
  flow_file = "./private-data/Yampa/2024uym_flow.csv",
  start_date = "2024-06-9",
  end_date = "2024-06-30",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)
load_uym24_event2$volume
load_uym24_event2$loads

load_uym24_event3 <- run_load_analysis(
  wq_file = "./private-data/Yampa/2024uym_wq.csv",
  flow_file = "./private-data/Yampa/2024uym_flow.csv",
  start_date = "2024-07-01",
  end_date = "2024-07-31",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)
load_uym24_event3$volume
load_uym24_event3$loads

total24_uym <- sum_load_objects(load_uym24_event1, load_uym24_event2, load_uym24_event3)
total24_uym$volume
total24_uym$loads
