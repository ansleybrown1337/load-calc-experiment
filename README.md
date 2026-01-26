# load-calc-experiment

A lightweight but robust R toolkit for calculating irrigation-event and seasonal water-quality loads by combining ISCO-style flow data with discrete grab or composite water-quality samples. The package is designed for real field datasets where logging intervals vary, analytes differ by unit, volumes may be missing, and paired-plot comparisons are required.

This repository grew out of edge-of-field irrigation studies (e.g., Fruita C, Upper Yampa, Gunnison, Molina) and is intended to be transparent, auditable, and easy to adapt to new sites.

---

## Key Capabilities

- Robust import of ISCO flow exports with variable headers and column counts
- Flexible datetime parsing with explicit time-zone handling
- Optional aggregation of flow to user-defined time steps
- Automatic volume estimation from average flow when volume is missing
- Water-quality preprocessing with:
  - analyte name normalization
  - non-detect handling
  - Fe and Se unit normalization (µg/L → mg/L)
- Event-based and seasonal load calculations
- Derived analytes (e.g., Total N = TKN + NO3-N + NO2-N)
- Explicit uncertainty bounds based on concentration standard deviation
- Summation of multiple irrigation events
- Conversion of loads to spatial units (lb/ac, kg/ha)
- Interactive Plotly diagnostics:
  - flow time series with bottle/sample markers
  - overlaid flow comparisons (C1 vs C2)
  - C1 vs C2 scatterplots on overlapping timesteps
- Automated C1 vs C2 comparison tables with non-overlap significance logic

---

## Repository Structure

```
load-calc-experiment/
│
├── code/
│   ├── load_exp.R              # Core functions
│   ├── fruita_C_2025.R         # Example paired-plot analysis
│   ├── upper-yampa-2022-2024.R # Multi-year watershed example
│   └── site-specific scripts
│
├── data/                        # Public example datasets
│
├── private-data/                # Ignored raw field data
│
├── figs/                        # Example outputs
│
├── README.md
├── LICENSE
└── .gitignore
```

`load_exp.R` is the only file required to use the toolkit. All other scripts are worked examples.

---

## Installation and Setup

Clone the repository and open the R project:

```r
source("./code/load_exp.R")
```

Required packages:

- dplyr
- tidyr
- lubridate
- plotly

All dependencies are loaded internally with informative error messages if missing.

---

## Core Workflow

### 1. Run a Load Analysis for a Single Event

```r
res <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date   = "2022-07-08",
  treatment_filter = "Outflow",
  user_interval = 4
)

res$volume   # total volume (L)
res$loads    # analyte loads (kg)
```

### 2. Combine Multiple Events

```r
total <- sum_load_objects(event1, event2, event3)
```

Volumes and analyte loads are summed explicitly. No implicit averaging occurs.

### 3. Convert to Spatial Units

```r
spatial <- convert_load_to_spatial(total, acreage = 8.2)

spatial$volume_acre_ft
spatial$loads
```

Outputs include lb/ac and kg/ha with upper and lower bounds.

---

## Paired-Plot and Treatment Comparisons

### Flow Diagnostics

Overlay two flow datasets in time:

```r
plot_avg_flow_timeseries_plotly(
  flow_file   = flowpath_C1,
  flow_file_2 = flowpath_C2,
  flow_units  = "cfs",
  flow_label_1 = "C1 (32-0-0 UAN)",
  flow_label_2 = "C2 (NPower 27-0-0)",
  plot_title  = "Fruita C: C1 vs C2 Flow"
)
```

### Flow Relationship on Overlapping Timesteps

```r
merged <- plot_c1_c2_flow_relationship_plotly(
  flow_file_C1 = flowpath_C1,
  flow_file_C2 = flowpath_C2,
  user_interval = 1,
  flow_units = "cfs"
)
```

This function:
- aligns datasets by datetime
- drops non-overlapping timesteps
- removes non-positive flows
- optionally adds 1:1 and linear fit lines

---

## C1 vs C2 Load Comparison Tables

```r
final_tbl <- make_c1_c2_comparison_table(
  c1_loads = spatial_loads_C1$loads,
  c2_loads = spatial_loads_C2$loads,
  c1_label = "C1: 32 UAN",
  c2_label = "C2: Npower"
)
```

The resulting table includes:

- total loads (lbs)
- spatial loads (lbs/ac)
- a qualitative significance flag

### Significance Logic

Significance is determined **only** from uncertainty bounds:

- `C1 > C2` if lower(C1) > upper(C2)
- `C2 > C1` if lower(C2) > upper(C1)
- `NS` otherwise

Point estimates are *not* used for significance decisions.

---

## Mathematical and Data Assumptions

1. **Concentration averaging**
   - All WQ samples within a window are averaged
   - Multiple sampling days trigger a warning but are allowed

2. **Non-detects**
   - Treated as zero by design (`<MDL → 0`)

3. **Uncertainty**
   - Bounds reflect ±1 SD of concentration only
   - If SD cannot be estimated (n < 2), bounds are NA

4. **Volumes**
   - Prefer measured volume from ISCO exports
   - Otherwise estimated from avg_flow × timestep
   - Optional volume override enables paired-plot normalization

5. **Fe and Se units**
   - Assumed reported in µg/L
   - Converted internally to mg/L before load calculation

6. **Derived analytes**
   - Total N = TKN + NO3-N + NO2-N
   - Uncertainty propagated via quadrature when available

7. **Time zones**
   - All datetimes parsed and stored in America/Denver by default

---

## Intended Use

This tool is designed for:

- edge-of-field irrigation studies
- paired treatment comparisons
- expert-witness and audit-ready analyses
- exploratory diagnostics prior to formal modeling

It is **not** a replacement for regulatory load estimation frameworks, but rather a transparent research-grade calculation engine.

---

## License

MIT License. Use freely with attribution.

---

## Author

Ansley J. Brown

Agricultural Water Quality Program
Colorado State University
