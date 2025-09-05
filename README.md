# Load Calculation Experiment

## Purpose

This script calculates water analyte loads based on **water quality (WQ)** and **flow** data, while also documenting important assumptions and limitations. It aggregates flow, maps analyte names to standardized abbreviations, and estimates uncertainty in loads using **standard deviations**. You can choose the treatment location (e.g., `Outflow`, `Inflow`) and set the flow aggregation interval.

## Author

**A.J. Brown**
*Agricultural Data Scientist, CSU Agricultural Water Quality Program*

## Features

* **Automated processing** of WQ and flow CSVs.
* **Custom treatment selection** via `treatment_filter`.
* **Analyte name mapping** to standardized abbreviations.
* **Flow aggregation** to user-defined intervals (e.g., 4-hour blocks).
* **Load calculation** (kg) with upper/lower bounds using ±1 SD (lower bound clamped at 0).
* **Clear status messages** during processing (import counts, date filtering, excluded analytes, etc.).
* **Simple outputs**: a list containing total volume (L) and a tidy analyte load table.

### What’s new (backend improvements)

* **EC25 & pH automatically excluded** from load results (with a status line).
* **Non-detects (ND and `<MDL`) treated as zeros** (commented as such in the code).
* **Prototype flow-weighted mean concentration (FWMC)** utilities included but **not active** (no behavior change unless you opt-in later).
* **Safer math**: volume sums use `na.rm=TRUE`; negative lower bounds clamped to 0; guards for zero/NA flow.

---

## Installation & Dependencies

Requires **R** and these libraries:

```r
install.packages(c("dplyr", "tidyr", "lubridate"))
```

---

## Usage

### 1) Import the functions

Save the script as `load_exp.R` and source it:

```r
source("./code/load_exp.R")
```

### 2) Run a load calculation

```r
load_results <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date = "2022-07-08",
  treatment_filter = "Outflow",  # or "Inflow", etc.
  user_interval = 4              # aggregation in hours
)
```

### 3) Access results

* **Total outflow volume (L)**

```r
load_results$volume
```

* **Analyte load table (kg)**

```r
load_results$loads
# Columns: analyte, load_kg, load_upper_kg, load_lower_kg
```

* **Access a specific analyte’s row**, e.g., NO3\_N:

```r
subset(load_results$loads, analyte == "NO3_N")
# or dplyr:
# dplyr::filter(load_results$loads, analyte == "NO3_N")
```

> Note: **EC25** and **pH** are always excluded from `load_results$loads` (status message is printed during the run).

---

## Units

* **Flow/Volume**
  Input volume (gallons) → converted to **liters (L)** using `1 gal = 3.78541 L`.
  Total and aggregated volumes are reported in **liters (L)**.

* **Concentrations**
  Expected to be **mg/L** (see WQ data format below).

* **Loads**
  Computed as `mg/L × L = mg`, then converted to **kg** by dividing by **1e6**.

* **Uncertainty bounds**
  `load_upper_kg` and `load_lower_kg` use ±1 SD around the mean concentration; lower bound is clamped to 0.

---

## Data Format Requirements

### Water Quality CSV

**Expected structure (headers on row 1):**

* `location.name`
* `treatment.name` (e.g., Outflow, Inflow, CT, MT, ST, W1, …)
* `event.type` (e.g., "Inflow", "Outflow", "Point Sample")
* `sample.id` (AWQP sample ID)
* `lab.id.x` (ALS Lab ID)
* `method` (e.g., "SM 4500-NH3", "SM 4500-NO3")
* `cas.number`
* `analyte` (e.g., "NITRATE AS N", "SELENIUM")
* `result` (**numeric concentration in mg/L**; see ND rules below)
* `units` (e.g., "MG/L", "UG/L")
* `collected` (MM/DD/YYYY HH\:MM)
* `received` (MM/DD/YYYY HH\:MM)

**ND / Censored values:**

* **By design, NDs are treated as zeros.**
* Strings like `"ND"`, `"nd"`, `"non-detect"`, `"n/d"`, and values like `"<0.01"` are parsed to **0**.
* If your lab provides MDLs separately and you prefer half-MDL or other rules, you can adjust the helper in `process_wq_data()`.

### Flow CSV

**Expected structure (headers start on row 7 in the source file; the code assigns names):**

* `datetime` (MM/DD/YYYY HH\:MM)
* `min_flow_gpm`
* `min_time` (HH\:MM\:SS AM/PM)
* `max_flow_gpm`
* `max_time` (HH\:MM\:SS AM/PM)
* `avg_flow_gpm`
* `volume_gal`
* `sample_event` (1 = Yes, NA = No)

> The script converts `volume_gal` → `volume_L` and aggregates by your specified `user_interval` (hours).

---

## Example Session Output (abridged)

During a run you’ll see messages like:

```
WQ Data Imported: 108 rows
Flow Data Imported: 39613 rows
WQ Data Processed
Flow Data Processed
Start Date: 2023-05-03 00:00:00 UTC
End Date: 2023-08-03 00:00:00 UTC
WQ Data Filtered: 24 rows
Flow Data Filtered: 10562 rows
Flow Data Aggregated: 552 rows
Total Volume, L: 1.11e+07 | Gallons: 2.94e+06
Note: The following analytes were excluded from load calculations because they are not mass-based: EC25, pH
Excluded non-mass analytes from load calc: EC25, pH
Final analytes in load table: NO2_N, NO3_N, PO4_P, TP, TDS, TKN, TSS, Se (g not kg), Fe (g not kg)
```

---

## Prototype: Flow-Weighted Mean Concentrations (FWMC)

The code includes **prototype** helpers to compute flow-weighted concentrations by aligning WQ timestamps to aggregated flow `time_block`s within a tolerance window.
**These functions are not active** in `run_load_analysis()` to avoid changing your current results.
To pilot later:

1. `wq_aligned <- align_wq_to_flow_prototype(wq, flow_aggregated)`
2. `analyte_summary <- compute_analyte_summary_fwmc_prototype(wq_aligned, flow_aggregated)`

---

## License

This project is licensed under the **GNU GPL v2**.

## Contact

For questions, improvements, or collaborations, contact **A.J. Brown** – CSU Agricultural Water Quality Program.
