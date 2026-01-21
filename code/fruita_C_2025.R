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

# Irrigation 5 sampled but not processed due to C1 power issues
# F25_irr5C1 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C1,
#   start_date = "2025-04-10 00:00",
#   end_date = "2025-04-11 10:30",
#   treatment_filter = "C1",  # User can now specify treatment!
#   user_interval = 1
# )

# Irrigation 6 not sampled, unsure why
# F25_irr6C1 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C1,
#   start_date = "2025-04-10 00:00",
#   end_date = "2025-04-11 10:30",
#   treatment_filter = "C1",  # User can now specify treatment!
#   user_interval = 1
# )

# ----- Irrigations - C2 -----

# Irrigation 1 - not captured in flow or wq data

# Irrigation 2
F25_irr2C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-05 00:00",
  end_date = "2025-05-07 00:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1
)

# Irrigation 3
F25_irr3C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-17 12:00",
  end_date = "2025-05-18 20:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1
)

# Irrigation 4
F25_irr4C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-29 10:00",
  end_date = "2025-05-30 20:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1
)

# Irrigation 5 sampled but not processed due to C1 power issues
# F25_irr5C2 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C2,
#   start_date = "2025-06-10 00:00",
#   end_date = "2025-06-11 12:00",
#   treatment_filter = "C2",  # User can now specify treatment!
#   user_interval = 1
# )

# Irrigation 6 not sampled, unsure why
# F25_irr6C2 <- run_load_analysis(
#   wq_file = wqpath,
#   flow_file = flowpath_C2,
#   start_date = "2025-06-21 00:00",
#   end_date = "2025-06-22 12:00",
#   treatment_filter = "C2",  # User can now specify treatment!
#   user_interval = 1
# )


# ----- Summary of Results -----

F25_irr1C1$volume
F25_irr1C1$loads
F25_irr2C1$volume
F25_irr2C1$loads
F25_irr3C1$volume
F25_irr3C1$loads
F25_irr4C1$volume
F25_irr4C1$loads

F25_irr2C2$volume
F25_irr2C2$loads
F25_irr3C2$volume
F25_irr3C2$loads
F25_irr4C2$volume
F25_irr4C2$loads



# C1: 32-0-0 UAN fertilizer
total_C1 <- sum_load_objects(#F25_irr1C1,
                             F25_irr2C1,
                             F25_irr3C1,
                             F25_irr4C1)
total_C1$volume
total_C1$loads

# C2: Npower 27-0-0 low vol
total_C2 <- sum_load_objects(F25_irr2C2,
                             F25_irr3C2,
                             F25_irr4C2)
total_C2$volume
total_C2$loads


# convert to spatial units
spatial_loads_C1 <- convert_load_to_spatial(total_C1, acreage = 8.2)
spatial_loads_C2 <- convert_load_to_spatial(total_C2, acreage = 9.9)


# View total volume in acre-feet
spatial_loads_C1$volume_acre_ft
spatial_loads_C2$volume_acre_ft


# View converted loads in lbs/acre
print(spatial_loads_C1$loads)
print(spatial_loads_C2$loads)


