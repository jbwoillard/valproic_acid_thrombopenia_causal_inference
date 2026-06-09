# Valproate Adverse Effects in Critically Ill Patients — Target Trial Emulation on MIMIC-IV

**Study type:** Observational target trial emulation  
**Data source:** MIMIC-IV v3.1 (Beth Israel Deaconess Medical Center)  
**Primary outcome:** Incident thrombocytopenia (platelets < 100 G/L, Days 0–30)  
**Exposure:** High-dose valproate (> 1000 mg/day vs. ≤ 1000 mg/day)  
**Methods:** G-formula, IPTW, AIPW (doubly robust), TMLE, DML, GRF, BART, BCF

---

## Data availability

**MIMIC-IV data are not included in this repository.**  
Access requires completion of a data use agreement with PhysioNet:  
<https://physionet.org/content/mimiciv/3.1/>

Once access is granted, download MIMIC-IV v3.1 and set the path in `R/config/config.R`:

```r
MIMIC_ROOT <- "/path/to/your/mimic-iv-3.1"
```

---

## Repository structure

```
.
├── R/
│   ├── config/
│   │   └── config.R                     # Paths, constants, study parameters
│   ├── utils/
│   │   └── utils.R                      # Logging, I/O helpers, SMD, winsorisation
│   ├── cohort/
│   │   ├── build_cohort.R               # Cohort construction (incident-user new-user design)
│   │   └── engineer_features.R          # Feature engineering, MICE imputation
│   ├── analysis/
│   │   ├── describe_population.R        # Table 1, descriptive statistics
│   │   ├── estimate_primary_ate.R       # G-formula, IPTW, AIPW (primary ATE)
│   │   ├── estimate_extended_ate.R      # TMLE, DML, BART, GRF (multi-method)
│   │   └── run_hte_grf.R                # Heterogeneous effects (GRF causal forest)
│   ├── sensitivity/
│   │   ├── run_vibration_analysis.R     # Vibration across analytical specifications
│   │   ├── run_secondary_outcomes.R     # Secondary outcomes (hepatotoxicity, hyperammonaemia)
│   │   ├── run_missing_data_sensitivity.R  # 4-scenario missing-data sensitivity
│   │   ├── run_grf_calibration.R        # GRF calibration tests (test_calibration)
│   │   ├── run_balance_diagnostics.R    # PS overlap, Love plots, ESS
│   │   ├── run_bcf_hte.R                # BCF Bayesian causal forest + AIPW cross-fitting
│   │   └── run_weight_normalized_sensitivity.R  # Pharmacological sensitivity (mg/kg/day)
│   └── concentration/
│       ├── describe_concentration_monitoring.R  # TDM coverage and selection bias
│       ├── build_concentration_subcohort.R      # Concentration-based subcohort
│       ├── estimate_concentration_ate.R         # Concentration-based ATE
│       ├── run_concentration_dose_response.R    # Continuous dose-response (splines)
│       ├── run_concentration_hte.R              # Concentration-based HTE
│       └── run_concentration_timing.R           # Trough timing reclassification
├── run_all.R                            # Master pipeline script
├── .gitignore
└── README.md
```

---

## How to run

### Prerequisites

Set the working directory to the project root before running any script:

```r
setwd("/path/to/project")
```

Or open the project in RStudio, which sets the working directory automatically.

### Step-by-step execution

Scripts must be run **in order** — each step depends on outputs from the previous one.

```r
source("run_all.R")  # Runs the full pipeline
```

Or run individual steps:

```r
source("R/cohort/build_cohort.R")
source("R/cohort/engineer_features.R")
source("R/analysis/describe_population.R")
source("R/analysis/estimate_primary_ate.R")
source("R/analysis/estimate_extended_ate.R")
source("R/analysis/run_hte_grf.R")
# ... sensitivity and concentration scripts as needed
```

---

## Expected outputs

All outputs are written to `outputs/` (not versioned):

```
outputs/
├── tables/        # CSV result tables
├── figures/       # PNG figures
├── diagnostics/   # Propensity score diagnostics
└── models/        # Saved causal forest objects
```

Intermediate data are written to `data_intermediate/` and final analytic datasets to `data_final/` (not versioned — contain derived MIMIC data).

---

## R package dependencies

```r
# Core
install.packages(c(
  "data.table", "arrow", "dplyr",
  "ggplot2", "patchwork", "knitr"
))

# Causal inference
install.packages(c(
  "WeightIt", "cobalt",    # IPTW / balance
  "grf",                   # GRF causal forest
  "lmtest", "sandwich"     # Robust inference
))

# Advanced estimators (optional for extended ATE)
install.packages(c(
  "tmle",       # TMLE
  "DoubleML",   # DML
  "dbarts",     # BART
  "bcf"         # Bayesian causal forest
))

# Imputation & missing data
install.packages("mice")

# Utilities
install.packages(c("here", "lubridate", "DBI", "duckdb"))

# Report rendering
install.packages(c("rmarkdown", "gtsummary", "splines"))
```

---

## Study design summary

| Component | Definition |
|-----------|------------|
| **Population** | Adult (≥ 18 y) incident VPA initiators during index hospitalisation |
| **Intervention** | Daily VPA dose > 1000 mg/day (high dose) |
| **Comparator** | Daily VPA dose ≤ 1000 mg/day (low dose) |
| **Time zero** | First VPA prescription (avoids immortal time bias) |
| **Follow-up** | 30 days or hospital discharge, whichever came first |
| **Primary outcome** | Incident thrombocytopenia (platelet count < 100 G/L) |
| **Secondary outcomes** | Thrombocytopenia < 50 G/L; ≥ 30% platelet decline; hepatotoxicity (ALT/AST > 3×ULN); hyperammonaemia (NH₃ > 55 µmol/L) |
| **Confounders** | 23 pre-time-zero covariates (demographics, severity, comorbidities, baseline labs, co-medications) |
| **Primary estimator** | AIPW (doubly robust; bootstrap 95% CI, 500 replications) |

---

## Notes on the concentration-based analysis

Scripts in `R/concentration/` implement an **exploratory** secondary analysis using measured plasma VPA concentrations (MIMIC-IV item 51008). This analysis is subject to severe positivity violations (IPTW ESS ≈ 3%) and selection bias from differential therapeutic drug monitoring. Results should not be interpreted as causal estimates.

---

## Citation

If you use this code, please cite the original article and the MIMIC-IV database:

> Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data*. 2023.

---

## License

Code is released under the MIT License. MIMIC-IV data are subject to their own data use agreement and are not redistributed here.
