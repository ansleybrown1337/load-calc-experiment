# working file to use the load_exp.R functions on real data

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2023 Gunnison Data

guni23_event1 <- run_load_analysis(
  wq_file = "./private-data/Gunnison/gunnison23_wq.csv",
  flow_file = "./private-data/Gunnison/gunnison23_flow.csv",
  start_date = "2023-05-03",
  end_date = "2023-08-03",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

guni23_event1$volume
guni23_event1$loads

guni23_event2 <- run_load_analysis(
  wq_file = "./private-data/Gunnison/gunnison23_wq.csv",
  flow_file = "./private-data/Gunnison/gunnison23_flow.csv",
  start_date = "2023-08-04",
  end_date = "2023-12-31",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

guni23_event2$volume
guni23_event2$loads


total23_guni <- sum_load_objects(guni23_event1, guni23_event2)
total23_guni$volume
total23_guni$loads

