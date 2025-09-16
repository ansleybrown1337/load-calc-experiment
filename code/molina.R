# working file to use the load_exp.R functions on real data

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2022 Molina Data - Data not reliable. Do not calculate!! - AJ 16 Sep 2025
# 
#   mol22_event1 <- run_load_analysis(
#     wq_file = "./private-data/Molina/molina22_wq_new.csv",
#     flow_file = "./private-data/Molina/molina22_flow.csv",
#     start_date = "2022-05-27",
#     end_date = "2022-06-07",
#     treatment_filter = "Outflow",  # User can now specify treatment!
#     user_interval = 4
#   )
#   
#   mol22_event1$volume
#   mol22_event1$loads
#   
#   mol22_event2 <- run_load_analysis(
#     wq_file = "./private-data/Molina/molina22_wq_new.csv",
#     flow_file = "./private-data/Molina/molina22_flow.csv",
#     start_date = "2022-05-27",
#     end_date = "2022-06-07",
#     treatment_filter = "Outflow",  # User can now specify treatment!
#     user_interval = 4
#   )
#   
#   mol22_event2$volume
#   mol22_event2$loads
#   
#   mol22_event3 <- run_load_analysis(
#     wq_file = "./private-data/Molina/molina22_wq_new.csv",
#     flow_file = "./private-data/Molina/molina22_flow.csv",
#     start_date = "2022-05-27",
#     end_date = "2022-06-07",
#     treatment_filter = "Outflow",  # User can now specify treatment!
#     user_interval = 4
#   )
#   
#   mol22_event3$volume
#   mol22_event3$loads
#   
#   total22_mol <- sum_load_objects(mol22_event1, mol22_event2, mol22_event3)
#   total22_mol$volume
#   total22_mol$loads
  
## 2023 Molina Data - bad flow data dates - Fixed 9/5/25 noted by AJ

mol23_event1 <- run_load_analysis(
  wq_file = "./private-data/Molina/molina23_wq.csv",
  flow_file = "./private-data/Molina/molina23_flow.csv",
  start_date = "2023-05-23",
  end_date = "2023-07-08",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

mol23_event1$volume
mol23_event1$loads

mol23_event2 <- run_load_analysis(
  wq_file = "./private-data/Molina/molina23_wq.csv",
  flow_file = "./private-data/Molina/molina23_flow.csv",
  start_date = "2023-07-09",
  end_date = "2023-09-30",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

mol23_event2$volume
mol23_event2$loads


total23_mol <- sum_load_objects(mol23_event1, mol23_event2)
total23_mol$volume
total23_mol$loads

## 2024 SC Data

mol24_event1 <- run_load_analysis(
  wq_file = "./private-data/Molina/molina24_wq.csv",
  flow_file = "./private-data/Molina/molina24_flow.csv",
  start_date = "2024-05-18",
  end_date = "2024-06-10",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

mol24_event1$volume
mol24_event1$loads

mol24_event2 <- run_load_analysis(
  wq_file = "./private-data/Molina/molina24_wq.csv",
  flow_file = "./private-data/Molina/molina24_flow.csv",
  start_date = "2024-06-16",
  end_date = "2024-07-13",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

mol24_event2$volume
mol24_event2$loads

total24_mol <- sum_load_objects(mol24_event1, mol24_event2)
total24_mol$volume
total24_mol$loads
