# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2025 UYM Data
flowpath <- "./private-data/UYM/2025/uym_flow_2025.csv"
wqpath <- "./private-data/UYM/2025/uym_wq_2025.csv"

# quick timeseries plot of flow data
plot_avg_flow_timeseries_plotly(flowpath,
                                flow_units = 'cfs')


# ----- Irrigation 1 -----
# not collected b/c of ISCO programming runtime error resulting in no sample

uym25_irr1 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-05-13 01:00", ## YYYY-MM-DD HH:MM
  end_date = "2025-05-22 12:00",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)


# ----- Irrigation 2 -----

uym25_irr2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-05-30 11:00", ## YYYY-MM-DD HH:MM
  end_date = "2025-06-03 00:00",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# ----- Irrigation 3 -----

uym25_irr3 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-10 13:00", ## YYYY-MM-DD HH:MM
  end_date = "2025-06-21 00:00",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# ----- Irrigation 4 -----

uym25_irr4 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath,
  start_date = "2025-06-29 16:00", ## YYYY-MM-DD HH:MM
  end_date = "2025-07-08 08:00",
  treatment_filter = "Outflow",  # User can now specify treatment!
  user_interval = 4
)

# ----- Summary of Results -----

uym25_irr1$volume
uym25_irr1$loads
uym25_irr2$volume
uym25_irr2$loads
uym25_irr3$volume
uym25_irr3$loads
uym25_irr4$volume
uym25_irr4$loads

# Total Loads for irr 2-4 (no data for irr 1)
total_uym <- sum_load_objects(uym25_irr2, uym25_irr3, uym25_irr4)
total_uym$volume
total_uym$loads

# convert to spatial units
spatial_loads <- convert_load_to_spatial(total_uym, acreage = 66)

# View total volume in acre-feet
spatial_loads$volume_acre_ft

# View converted loads in lbs/acre
print(spatial_loads$loads)


