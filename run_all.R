# =============================================================================
# run_all.R — Full analysis pipeline
#
# Requirements:
#   1. MIMIC-IV v3.1 data accessible at a path configured in R/config/config.R
#   2. R packages installed (see README.md for the full list)
#   3. Working directory set to the project root (folder containing run_all.R)
#
# Usage:
#   Rscript run_all.R
#   or interactively: source("run_all.R")
#
# Each step can also be sourced individually. Steps must be run in order
# since each step depends on outputs from the previous one.
# =============================================================================

source("R/config/config.R")
source("R/utils/utils.R")

run_step <- function(script, label) {
  cat(sprintf("\n[STEP] %s\n", label))
  t0 <- proc.time()
  tryCatch(
    source(script, local = FALSE),
    error = function(e) cat(sprintf("[ERROR] %s failed: %s\n", label, e$message))
  )
  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("[DONE] %s | %.0f sec\n", label, elapsed))
}

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("VPA MIMIC-IV — Target Trial Emulation Pipeline\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# ---- Core pipeline --------------------------------------------------------

run_step("R/cohort/build_cohort.R",          "1 | Build analytic cohort")
run_step("R/cohort/engineer_features.R",     "2 | Engineer features & MICE imputation")

# ---- Descriptive ----------------------------------------------------------

run_step("R/analysis/describe_population.R", "3 | Descriptive statistics")

# ---- Primary causal analysis ---------------------------------------------

run_step("R/analysis/estimate_primary_ate.R",   "4 | Primary ATE (G-formula, IPTW, AIPW)")
run_step("R/analysis/estimate_extended_ate.R",  "5 | Extended ATE (TMLE, DML, BART, GRF)")
run_step("R/analysis/run_hte_grf.R",            "6 | Heterogeneous treatment effects (GRF)")

# ---- Sensitivity analyses ------------------------------------------------

run_step("R/sensitivity/run_vibration_analysis.R",           "7  | Vibration analysis")
run_step("R/sensitivity/run_secondary_outcomes.R",           "8  | Secondary outcomes")
run_step("R/sensitivity/run_missing_data_sensitivity.R",     "9  | Missing-data sensitivity")
run_step("R/sensitivity/run_grf_calibration.R",              "10 | GRF calibration")
run_step("R/sensitivity/run_balance_diagnostics.R",          "11 | Balance & overlap diagnostics")
run_step("R/sensitivity/run_bcf_hte.R",                      "12 | BCF HTE + AIPW cross-fitting")
run_step("R/sensitivity/run_weight_normalized_sensitivity.R","13 | Weight-normalised dose sensitivity")

# ---- Concentration-based analysis (exploratory) --------------------------

cat("\n--- Concentration-based exploratory analyses ---\n")
run_step("R/concentration/describe_concentration_monitoring.R", "C1 | TDM monitoring characterisation")
run_step("R/concentration/build_concentration_subcohort.R",     "C2 | Concentration subcohort")
run_step("R/concentration/estimate_concentration_ate.R",        "C3 | Concentration-based ATE")
run_step("R/concentration/run_concentration_dose_response.R",   "C4 | Dose-response splines")
run_step("R/concentration/run_concentration_hte.R",             "C5 | Concentration-based HTE")
run_step("R/concentration/run_concentration_timing.R",          "C6 | Trough timing analysis")


