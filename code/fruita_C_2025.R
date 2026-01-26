# working file to use the load_exp.R functions on real data

# set working directory
setwd("C:/Users/AJ-CPU/Documents/GitHub/load-calc-experiment")

# Load the functions from load_exp.R
source("./code/load_exp.R")

## 2025 Fruita C1 Data
flowpath_C1 <- "./private-data/Fruita C/fruitaC1_flow.csv"
flowpath_C2 <- "./private-data/Fruita C/fruitaC2_flowfull.csv"
wqpath <- "./private-data/Fruita c/fruitaC_wq.csv"

# quick timeseries plot of flow data
plot_avg_flow_timeseries_plotly(flowpath_C1,
                                flow_units = 'cfs',
                                plot_title = "Fruita C1 Flow Timeseries")
plot_avg_flow_timeseries_plotly(flowpath_C2,
                                flow_units = 'cfs',
                                plot_title = "Fruita C2 Flow Timeseries")

# C1 vs C2 comparison (overlay)
plot_avg_flow_timeseries_plotly(
  flow_file   = flowpath_C1,
  flow_file_2 = flowpath_C2,
  flow_units  = "cfs",
  flow_label_1 = "C1 (32-0-0 UAN)",
  flow_label_2 = "C2 (NPower 27-0-0)",
  plot_title  = "Fruita C: C1 vs C2 Flow Timeseries"
)

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

# Irrigation 1 - not captured in flow or wq data; C2 ramp flume blew out

# Irrigation 2
F25_irr2C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-05 00:00",
  end_date = "2025-05-07 00:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1,
  volume_override_L = F25_irr2C1$volume  # use C1 volume for C2 load calc
)

# Irrigation 3
F25_irr3C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-17 12:00",
  end_date = "2025-05-18 20:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1,
  volume_override_L = F25_irr3C1$volume  # use C1 volume for C2 load calc
)

# Irrigation 4
F25_irr4C2 <- run_load_analysis(
  wq_file = wqpath,
  flow_file = flowpath_C2,
  start_date = "2025-05-29 10:00",
  end_date = "2025-05-30 20:00",
  treatment_filter = "C2",  # User can now specify treatment!
  user_interval = 1,
  volume_override_L = F25_irr4C1$volume  # use C1 volume for C2 load calc
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

# sig diff fxn for comparing between load objects
final_tbl <- make_c1_c2_comparison_table(
  c1_loads = spatial_loads_C1$loads,
  c2_loads = spatial_loads_C2$loads,
  c1_label = "C1: 32 UAN",
  c2_label = "C2: Npower"
)

# final_tbl is a tibble you can also write to CSV if you want:
# write.csv(final_tbl, "c1_c2_load_comparison.csv", row.names = FALSE)



# View converted loads in lbs/acre
print(spatial_loads_C1$loads)
print(spatial_loads_C2$loads)


# ------ c1 v c2 plotter -------
plot_c1_c2_flow_relationship_plotly <- function(flow_file_C1,
                                                flow_file_C2,
                                                start_date = NULL,
                                                end_date   = NULL,
                                                tz = "America/Denver",
                                                user_interval = NULL,     # NULL = no aggregation; else hours (e.g., 1)
                                                flow_units = c("gpm", "cfs", "L/s"),
                                                plot_title = "C1 vs C2 flow relationship",
                                                include_fit = TRUE,
                                                include_1to1 = TRUE) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install it with install.packages('plotly').", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  
  flow_units <- match.arg(flow_units)
  
  # ---- unit conversion (from gpm) ----
  # 1 cfs = 448.831 gpm
  # 1 gpm = 3.78541 L/min = 0.0630901667 L/s
  conv_factor <- switch(
    flow_units,
    "gpm" = 1,
    "cfs" = 1 / 448.831,
    "L/s" = 3.78541 / 60
  )
  unit_label <- switch(flow_units, "gpm" = "gpm", "cfs" = "cfs", "L/s" = "L/s")
  
  # Parse optional window
  start_dt <- if (!is.null(start_date)) as.POSIXct(trimws(start_date), tz = tz) else NULL
  end_dt   <- if (!is.null(end_date))   as.POSIXct(trimws(end_date),   tz = tz) else NULL
  
  # ---- Load + process ----
  c1_raw <- import_flow_data_flex(flow_file_C1, tz = tz)
  c2_raw <- import_flow_data_flex(flow_file_C2, tz = tz)
  
  c1 <- process_flow_data_flex(c1_raw, tz = tz)
  c2 <- process_flow_data_flex(c2_raw, tz = tz)
  
  # ---- Optional window filtering ----
  if (!is.null(start_dt)) {
    c1 <- dplyr::filter(c1, datetime >= start_dt)
    c2 <- dplyr::filter(c2, datetime >= start_dt)
  }
  if (!is.null(end_dt)) {
    c1 <- dplyr::filter(c1, datetime <= end_dt)
    c2 <- dplyr::filter(c2, datetime <= end_dt)
  }
  
  if (nrow(c1) == 0 || nrow(c2) == 0) {
    warning("After filtering, one of the datasets has 0 rows. Nothing to compare.", call. = FALSE)
    return(invisible(NULL))
  }
  
  # ---- Optional aggregation to common interval ----
  # Uses your existing aggregate_flow_data() which returns time_block, avg_flow_gpm, volume_L, etc.
  if (!is.null(user_interval) && is.finite(user_interval) && user_interval > 0) {
    c1a <- aggregate_flow_data(c1, user_interval = user_interval) %>%
      dplyr::rename(datetime = time_block, avg_flow_gpm_C1 = avg_flow_gpm)
    c2a <- aggregate_flow_data(c2, user_interval = user_interval) %>%
      dplyr::rename(datetime = time_block, avg_flow_gpm_C2 = avg_flow_gpm)
  } else {
    c1a <- c1 %>%
      dplyr::select(datetime, avg_flow_gpm) %>%
      dplyr::rename(avg_flow_gpm_C1 = avg_flow_gpm)
    c2a <- c2 %>%
      dplyr::select(datetime, avg_flow_gpm) %>%
      dplyr::rename(avg_flow_gpm_C2 = avg_flow_gpm)
  }
  
  # ---- Merge on datetime; drop non-overlap automatically (inner join) ----
  merged <- dplyr::inner_join(c1a, c2a, by = "datetime") %>%
    dplyr::filter(
      is.finite(avg_flow_gpm_C1),
      is.finite(avg_flow_gpm_C2),
      avg_flow_gpm_C1 > 0,
      avg_flow_gpm_C2 > 0
    ) %>%
    dplyr::mutate(
      flow_C1 = avg_flow_gpm_C1 * conv_factor,
      flow_C2 = avg_flow_gpm_C2 * conv_factor
    )
  
  if (nrow(merged) == 0) {
    warning(
      paste0(
        "No overlapping datetimes between C1 and C2 after filtering/aggregation.\n",
        "Tip: set user_interval (e.g., 1) to align irregular logging intervals."
      ),
      call. = FALSE
    )
    return(invisible(NULL))
  }
  
  # ---- Build plot ----
  hover_text <- paste0(
    "Datetime: ", merged$datetime,
    "<br>C1: ", sprintf("%.3f", merged$flow_C1), " ", unit_label,
    "<br>C2: ", sprintf("%.3f", merged$flow_C2), " ", unit_label
  )
  
  p <- plotly::plot_ly(
    data = merged,
    x = ~ flow_C1,
    y = ~ flow_C2,
    type = "scatter",
    mode = "markers",
    text = hover_text,
    hoverinfo = "text",
    name = "Overlapping timesteps"
  ) %>%
    plotly::layout(
      title = plot_title,
      xaxis = list(title = paste0("C1 flow (", unit_label, ")")),
      yaxis = list(title = paste0("C2 flow (", unit_label, ")")),
      hovermode = "closest"
    )
  
  # ---- Optional 1:1 line ----
  if (isTRUE(include_1to1)) {
    rng <- range(c(merged$flow_C1, merged$flow_C2), na.rm = TRUE)
    line_df <- data.frame(x = rng, y = rng)
    
    p <- p %>% plotly::add_trace(
      data = line_df,
      x = ~ x, y = ~ y,
      type = "scatter",
      mode = "lines",
      name = "1:1 line",
      inherit = FALSE,
      hoverinfo = "skip"
    )
  }
  
  # ---- Optional linear fit line (C2 ~ C1) ----
  if (isTRUE(include_fit) && nrow(merged) >= 3) {
    fit <- stats::lm(flow_C2 ~ flow_C1, data = merged)
    co <- stats::coef(fit)
    
    x_fit <- seq(min(merged$flow_C1, na.rm = TRUE), max(merged$flow_C1, na.rm = TRUE), length.out = 100)
    y_fit <- co[[1]] + co[[2]] * x_fit
    fit_df <- data.frame(x = x_fit, y = y_fit)
    
    fit_label <- paste0("Fit: C2 = ", signif(co[[1]], 4), " + ", signif(co[[2]], 4), "·C1")
    
    p <- p %>% plotly::add_trace(
      data = fit_df,
      x = ~ x, y = ~ y,
      type = "scatter",
      mode = "lines",
      name = fit_label,
      inherit = FALSE,
      hoverinfo = "skip"
    )
  }
  
  print(p)
  
  attr(merged, "flow_units") <- flow_units
  invisible(merged)
}

# scatterplot of C1 vs C2 flow to find relationship
# C1 vs C2 relationship during shared irrigations (example window)
merged_cc <- plot_c1_c2_flow_relationship_plotly(
  flow_file_C1 = flowpath_C1,
  flow_file_C2 = flowpath_C2,
  start_date = "2025-04-05 00:00",
  end_date   = "2025-06-30 23:59",
  tz = "America/Denver",
  user_interval = 1,       # strongly recommended for alignment
  flow_units = "cfs",
  plot_title = "Fruita C: C1 vs C2 flow (overlapping timesteps)"
)

# If you want the merged data for export or modeling:
head(merged_cc)


