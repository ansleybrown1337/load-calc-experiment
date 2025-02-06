# working file to use the load_exp.R functions on real data

# Load the functions from load_exp.R
source("./code/load_exp.R")

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

load_uym23_event1 <- run_load_analysis(
  wq_file = "./private-data/Yampa 2023/2023uym_wq.csv",
  flow_file = "./private-data/Yampa 2023/2023uym_flow.csv",
  start_date = "2023-06-01",
  end_date = "2022-06-14",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# Access results
load_uym23_event1$volume
load_uym23_event1$loads
