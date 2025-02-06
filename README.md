# Load Calculation Experiment

## Purpose
This script calculates water analyte loads based on **water quality data (WQ)** and **flow data**, while also showcasing the limitations of the methods used. It aggregates flow rates, maps analyte names to standardized abbreviations, and estimates uncertainty in load calculations using **standard deviations**. The user can specify the treatment location (e.g., `Outflow`, `Inflow`) and set custom aggregation intervals for flow data.

## Author
**A.J. Brown**  
_Agricultural Data Scientist, CSU Agricultural Water Quality Program_

## Features
- **Automated Processing:** Imports and processes WQ and flow data.
- **Custom Treatment Selection:** Users can specify `treatment.name` dynamically.
- **Analyte Mapping:** Matches old and new analyte names to standardized abbreviations.
- **Flow Aggregation:** Supports user-defined time intervals (e.g., 4-hour blocks).
- **Load Calculation:** Computes analyte loads in kg with upper/lower bounds.
- **Easy-to-Access Output:** Returns a structured list containing volume and analyte loads.

---

## Installation & Dependencies
This script requires **R** and the following libraries:
```r
install.packages(c("dplyr", "tidyr"))
```

## Usage
### 1. Import the Script
Save the script as `load_calc.R` and source it in R:
```r
source("load_calc.R")
```

### 2. Run Load Calculation
To execute the script, use:
```r
load_results <- run_load_analysis(
  wq_file = "./data/2022uym_wq.csv",
  flow_file = "./data/2022uym_flow.csv",
  start_date = "2022-06-12",
  end_date = "2022-07-08",
  treatment_filter = "Outflow",  # Change to "Inflow" or other treatments as needed
  user_interval = 4  # Time aggregation in hours
)
```

### 3. Access Results
#### a) Total Water Volume (L)
```r
load_results$volume
```
#### b) Analyte Load Estimates
```r
load_results$loads
```
#### c) Accessing a Specific Analyte Load
```r
load_results$NO3_N$load_kg        # Mean estimated load
load_results$NO3_N$load_upper_kg  # Upper bound
load_results$NO3_N$load_lower_kg  # Lower bound
```

## Units

## Data Format Requirements
### **Water Quality Data (CSV)**
**Expected structure before processing:**
- **Headers start on row 1**
- **Columns:**
  - `location.name` (e.g. "Kerbel")
  - `treatment.name` (e.g., Outflow, Inflow, CT, MT, ST, W1, etc.)
  - `event.type` (e.g., "Inflow", "Outflow", "Point Sample")
  - `sample.id` (i.e., AWQP sample ID)
  - `lab.id.x` (i.e., ALS Lab ID)
  - `method` (e.g., "SM 4500-NH3", "SM 4500-NO3")
  - `cas.number` (i.e., ALS Lab Information)
  - `analyte` (e.g., "NITRATE AS N", "SELENIUM")
  - `result` (numeric concentration value)
  - `units` (e.g., "MG/L", "UG/L")
  - `collected` (datetime, format: MM/DD/YYYY HH:MM)
  - `received` (datetime, format: MM/DD/YYYY HH:MM)

### **Flow Data (CSV)**
**Expected structure before processing:**
- **Headers start on row 7**
- **Columns (not named in file due to no headers, but named in this code):**
  - `datetime` (timestamp of flow measurement, format: MM/DD/YYYY HH:MM)
  - `min_flow_gpm` (minimum flow in gpm)
  - `min_time` (time when min flow occurred, format: HH:MM:SS AM/PM)
  - `max_flow_gpm` (maximum flow in gpm)
  - `max_time` (time when max flow occurred, format: HH:MM:SS AM/PM)
  - `avg_flow_gpm` (average flow in gpm)
  - `volume_gal` (total volume in gallons)
  - `sample_event` (indicator if a sample was collected, 1 = Yes, NA = No)

### Output Units in R Object

- Total volume is converted to liters (L) (1 gallon = 3.78541 L).
- Aggregated volumes are reported in liters (L).
- Final loads are converted to kilograms (kg) (mg / 1e6).

---

## Example R Code
```r
> load_results <- run_load_analysis(
+   wq_file = "./data/2022uym_wq.csv",
+   flow_file = "./data/2022uym_flow.csv",
+   start_date = "2022-06-12",
+   end_date = "2022-07-08",
+   treatment_filter = "Outflow",  # User can now specify treatment!
+   user_interval = 4
+ )
> # Access results
> load_results$volume
[1] 70124942
> load_results$loads
  analyte     load_kg load_upper_kg load_lower_kg
1    EC25         NaN           NaN           NaN
2   NO2_N   91.162425     91.162425    91.1624252
3   NO3_N    0.000000      0.000000     0.0000000
4      pH         NaN           NaN           NaN
5   PO4_P    0.000000      0.000000     0.0000000
6      Se    0.000000      0.000000     0.0000000
7      TP    1.928436      4.155938    -0.2990661
8     TSS 1869.998465   2103.748274  1636.2486571
> load_results$NO3_N$load_kg
[1] 0
> load_results$NO3_N$load_upper_kg
[1] 0
> load_results$NO3_N$load_lower_kg
[1] 0
```

---

## License
This project is licensed under the **GNU GPL v2 License**.

## Contact
For questions, improvements, or collaborations, contact **A.J. Brown** at **CSU Agricultural Water Quality Program**.

