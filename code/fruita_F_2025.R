# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2025 Fruita F Data
flowpath <- "./private-data/Fruita F/fruitaF_flow.csv"
wqpath <- "./private-data/Fruita F/fruitaF_wq_2025.csv"

# quick timeseries plot of flow data
plot_avg_flow_timeseries_plotly(flowpath)


# ----- Irrigation 1 - Sets 1-4 -----

F25_irr1set1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-04-10 00:00",
  end_date = "2025-04-11 10:30",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr1set2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-04-11 10:45",
  end_date = "2025-04-12 10:30",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr1set3 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-04-12 10:45",
  end_date = "2025-04-13 10:30 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr1set4 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-04-13 10:45",
  end_date = "2025-04-14 10:30 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

# ----- Irrigation 2 - not collected -----

# ----- Irrigation 3 - Sets 1-4 -----

F25_irr3set1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-09 16:00",
  end_date = "2025-06-10 10:30 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr3set2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-10 13:30",
  end_date = "2025-06-11 09:00 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr3set3 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-11 09:15",
  end_date = "2025-06-12 09:00 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

F25_irr3set4 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-12 09:15",
  end_date = "2025-06-13 09:30 ",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 1
)

# ----- Summary of Results -----

F25_irr1set1$volume
F25_irr1set1$loads
F25_irr1set2$volume
F25_irr1set2$loads
F25_irr1set3$volume
F25_irr1set3$loads
F25_irr1set4$volume
F25_irr1set4$loads

F25_irr3set1$volume
F25_irr3set1$loads
F25_irr3set2$volume
F25_irr3set2$loads
F25_irr3set3$volume
F25_irr3set3$loads
F25_irr3set4$volume
F25_irr3set4$loads

# F1: 28-0-0 fertilizer
total_F1 <- sum_load_objects(F25_irr1set1, F25_irr3set1)
total_F1$volume
total_F1$loads

# F2: 28-0-0 fertilizer (same as F1)
total_F2 <- sum_load_objects(F25_irr1set2, F25_irr3set2)
total_F2$volume
total_F2$loads

# F3: Triple N fertilizer
total_F3 <- sum_load_objects(F25_irr1set3, F25_irr3set3)
total_F3$volume
total_F3$loads

# F4: 32-0-0 UAN fertilizer
total_F4 <- sum_load_objects(F25_irr1set4, F25_irr3set4)
total_F4$volume
total_F4$loads

# convert to spatial units
# Assume we have a load object and field size of 50 acres
spatial_loads_F1 <- convert_load_to_spatial(total_F1, acreage = 13.2)
spatial_loads_F2 <- convert_load_to_spatial(total_F2, acreage = 12.9)
spatial_loads_F3 <- convert_load_to_spatial(total_F3, acreage = 13.1)
spatial_loads_F4 <- convert_load_to_spatial(total_F4, acreage = 11.6)

# View total volume in acre-feet
spatial_loads_F1$volume_acre_ft
spatial_loads_F2$volume_acre_ft
spatial_loads_F3$volume_acre_ft
spatial_loads_F4$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads_F1$loads)
print(spatial_loads_F2$loads)
print(spatial_loads_F3$loads)
print(spatial_loads_F4$loads)

