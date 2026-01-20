# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2025 Fruita C1 Data
flowpath_C1 <- "./private-data/Fruita C/fruitaC1_flow.csv"
flowpath_C2 <- "./private-data/Fruita C/fruitaC2_flow.csv"
wqpath <- "./private-data/Fruita c/fruitaC_wq.csv"

# quick timeseries plot of flow data
plot_avg_flow_timeseries_plotly(flowpath_C1,
                                flow_units = 'cfs')
plot_avg_flow_timeseries_plotly(flowpath_C2,
                                flow_units = 'cfs')


# ----- Irrigations - C1 -----
# Set 3 of each irrigation was the plot before the treatment and used here
# for comparison against C2

# Irrigation 1 erroneously sampled first and set, not 3rd set (as intended)
# b/c this second set is still control plot (32-0-0) we will still use it
F25_irr1C1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C1,
  start_date = "2025-04-14 22:00",
  end_date = "2025-04-16 00:00",
  treatment_filter = "C1",  # User can now specify treatment!
  user_interval = 1
)

F25_irr2C1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C1,
  start_date = "2025-05-05 00:00",
  end_date = "2025-05-06 00:00",
  treatment_filter = "C1",  # User can now specify treatment!
  user_interval = 1
)

F25_irr3C1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C1,
  start_date = "2025-05-16 22:00",
  end_date = "2025-05-17 14:00",
  treatment_filter = "C1",  # User can now specify treatment!
  user_interval = 1
)

F25_irr4C1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C1,
  start_date = "2025-05-28 21:00",
  end_date = "2025-05-29 13:00",
  treatment_filter = "C1",  # User can now specify treatment!
  user_interval = 1
)

# Irrigation 5 not processed due to power issues
# F25_irr5C1 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C1,
#   start_date = "2025-04-10 00:00",
#   end_date = "2025-04-11 10:30",
#   treatment_filter = "C1",  # User can now specify treatment!
#   user_interval = 1
# )

# not processed, unsure why
# F25_irr6C1 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C1,
#   start_date = "2025-04-10 00:00",
#   end_date = "2025-04-11 10:30",
#   treatment_filter = "C1",  # User can now specify treatment!
#   user_interval = 1
# )

# ----- Irrigations - C2 -----

F25_irr1C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-04-10 00:00",
  end_date = "2025-04-11 10:30",
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

