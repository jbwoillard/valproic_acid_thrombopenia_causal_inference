# =============================================================================
# run_secondary_outcomes.R
# Step 8 — Secondary outcomes and raw comparisons
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table); library(arrow); library(grf)

log_msg("=== STARTING ===")
t0 <- timer_start()

ensure_dirs(TABLES_DIR)

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
ate_results <- fread(file.path(TABLES_DIR, "ate_results.csv"))
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
cov_vars <- intersect(propensity_vars, names(cohort))

cohort_clean <- copy(cohort)
for (v in cov_vars) {
  if (is.numeric(cohort_clean[[v]])) {
    cohort_clean[is.na(get(v)), (v) := median(cohort_clean[[v]], na.rm=TRUE)]
  }
}

# =============================================================================
# A4 — Comparaisons brutes avec tests (chi-carré / Fisher)
# =============================================================================
log_msg("A4. Comparaisons brutes des outcomes...")

run_crude_test <- function(outcome_col, label, cohort_dt) {
  Y1 <- as.integer(cohort_dt[[outcome_col]][cohort_dt$high_dose==1])
  Y0 <- as.integer(cohort_dt[[outcome_col]][cohort_dt$high_dose==0])
  Y1 <- Y1[!is.na(Y1)]; Y0 <- Y0[!is.na(Y0)]
  n1 <- length(Y1); n0 <- length(Y0)
  ev1 <- sum(Y1); ev0 <- sum(Y0)
  p1 <- ev1/n1; p0 <- ev0/n0
  rd <- p1 - p0
  rr <- if(p0 > 0) p1/p0 else NA_real_

  # Test chi-carré (ou Fisher si petits effectifs)
  ct <- tryCatch({
    mat <- matrix(c(ev1, n1-ev1, ev0, n0-ev0), nrow=2)
    if (min(mat) < 5) fisher.test(mat) else chisq.test(mat, correct=FALSE)
  }, error=function(e) list(p.value=NA_real_))

  data.table(
    outcome       = label,
    n_high        = n1,
    n_low         = n0,
    events_high   = ev1,
    events_low    = ev0,
    pct_high      = round(p1*100, 1),
    pct_low       = round(p0*100, 1),
    rd_crude      = round(rd*100, 2),
    rr_crude      = round(rr, 2),
    p_value       = round(ct$p.value, 4),
    test_method   = if(min(ev1, ev0, n1-ev1, n0-ev0) < 5) "Fisher" else "Chi-carré"
  )
}

outcomes_to_test <- list(
  list(col="outcome_thrombo_100", label="Thrombopénie < 100 G/L"),
  list(col="outcome_thrombo_50",  label="Thrombopénie < 50 G/L"),
  list(col="outcome_thrombo_rel30", label="Baisse plaquettes > 30%"),
  list(col="outcome_hepato_3xULN", label="Hépatotoxicité > 3xULN"),
  list(col="outcome_hypernh3",    label="Hyperammoniémie > 55 µmol/L"),
  list(col="outcome_composite",   label="Outcome composite")
)

raw_comparisons <- rbindlist(lapply(outcomes_to_test, function(oc) {
  tryCatch(run_crude_test(oc$col, oc$label, cohort_clean),
           error=function(e) data.table(outcome=oc$label, n_high=NA, n_low=NA,
                                         events_high=NA, events_low=NA, pct_high=NA, pct_low=NA,
                                         rd_crude=NA, rr_crude=NA, p_value=NA, test_method=NA))
}))

safe_write_csv(raw_comparisons, file.path(TABLES_DIR, "raw_outcome_comparisons.csv"), overwrite=TRUE)
log_msg("Comparaisons brutes:")
print(raw_comparisons[, .(outcome, pct_high, pct_low, rd_crude, p_value, test_method)])

# =============================================================================
# A5 — Intégration ATE GRF dans tableau principal
# =============================================================================
log_msg("A5. Intégration ATE GRF dans tableau principal...")

grf_ate <- tryCatch(fread(file.path(TABLES_DIR, "grf_ate.csv")), error=function(e) NULL)

# Construire tableau ATE étendu avec GRF
if (!is.null(grf_ate)) {
  ate_grf_row <- data.table(
    estimator   = "7. GRF Causal Forest (ATE global)",
    outcome     = "Thrombopénie < 100 G/L (J0-J30)",
    estimate_rd = round(grf_ate$estimate, 4),
    ci_lower    = round(grf_ate$ci_lower, 4),
    ci_upper    = round(grf_ate$ci_upper, 4),
    note        = sprintf("ATE GRF — IC Neyman, 2000 arbres | Calibration: diff_p=0.894 (HTE non démontrée)")
  )
  ate_updated <- rbind(
    ate_results[outcome == "Thrombopénie < 100 G/L (J0-J30)"],
    ate_grf_row,
    fill=TRUE
  )
  log_msg(sprintf("ATE GRF intégré: RD=%.4f [%.4f, %.4f]",
                  grf_ate$estimate, grf_ate$ci_lower, grf_ate$ci_upper))
} else {
  ate_updated <- ate_results
  log_msg("WARN: grf_ate.csv non disponible")
}

# =============================================================================
# A11 — Outcomes secondaires dose-based avec IC
# =============================================================================
log_msg("A11. Outcomes secondaires dose-based (AIPW + bootstrap 200 reps)...")

ps_formula_d <- as.formula(paste("high_dose ~", paste(cov_vars, collapse="+")))
ps_m_d <- glm(ps_formula_d, data=as.data.frame(cohort_clean), family=binomial)
ps_hat_d <- pmax(pmin(predict(ps_m_d, type="response"), 0.99), 0.01)
A <- cohort_clean$high_dose
n <- nrow(cohort_clean)

run_aipw_with_ci <- function(outcome_col, label, B=200) {
  Y <- as.integer(cohort_clean[[outcome_col]])
  if (sum(!is.na(Y)) < 30 || sum(Y, na.rm=TRUE) < 5) {
    return(data.table(estimator="AIPW", outcome=label, estimate_rd=NA,
                      ci_lower=NA, ci_upper=NA, note="Trop peu de cas"))
  }
  f_oc <- as.formula(paste(outcome_col, "~ high_dose +", paste(cov_vars, collapse="+")))
  m_oc <- tryCatch(glm(f_oc, data=as.data.frame(cohort_clean), family=binomial),
                    error=function(e) NULL)
  if (is.null(m_oc)) return(data.table(estimator="AIPW", outcome=label, estimate_rd=NA,
                                        ci_lower=NA, ci_upper=NA, note="Erreur modèle"))
  d1 <- copy(cohort_clean); d1[, high_dose := 1L]
  d0 <- copy(cohort_clean); d0[, high_dose := 0L]
  mu1 <- predict(m_oc, newdata=d1, type="response")
  mu0 <- predict(m_oc, newdata=d0, type="response")
  eif1 <- mu1 + A/ps_hat_d * (Y - mu1)
  eif0 <- mu0 + (1-A)/(1-ps_hat_d) * (Y - mu0)
  ate_pt <- mean(eif1) - mean(eif0)

  # Bootstrap léger
  set.seed(RANDOM_SEED + 99)
  boot_v <- sapply(seq_len(B), function(b) {
    tryCatch({
      idx <- sample(n, n, replace=TRUE)
      bd <- cohort_clean[idx,]
      ps_b <- pmax(pmin(predict(glm(ps_formula_d, data=as.data.frame(bd), family=binomial),
                                 type="response"), 0.99), 0.01)
      m_b <- glm(f_oc, data=as.data.frame(bd), family=binomial)
      b1 <- copy(bd); b1[, high_dose := 1L]
      b0 <- copy(bd); b0[, high_dose := 0L]
      m1b <- predict(m_b, newdata=b1, type="response")
      m0b <- predict(m_b, newdata=b0, type="response")
      Ab <- bd$high_dose; Yb <- as.integer(bd[[outcome_col]])
      mean(m1b + Ab/ps_b*(Yb-m1b)) - mean(m0b + (1-Ab)/(1-ps_b)*(Yb-m0b))
    }, error=function(e) NA_real_)
  })
  ci <- quantile(boot_v, c(0.025, 0.975), na.rm=TRUE)
  n_ev <- sum(Y, na.rm=TRUE)
  log_msg(sprintf("  %s: N events=%d, ATE=%.4f [%.4f, %.4f]", label, n_ev, ate_pt, ci[1], ci[2]))
  data.table(estimator="AIPW", outcome=label,
             estimate_rd=round(ate_pt, 4),
             ci_lower=round(ci[1], 4),
             ci_upper=round(ci[2], 4),
             note=sprintf("IC bootstrap %d reps", B))
}

sec_outcomes <- list(
  list(col="outcome_hepato_3xULN", label="Hépatotoxicité ALT/AST > 3xULN"),
  list(col="outcome_hypernh3",     label="Hyperammoniémie NH3 > 55 µmol/L"),
  list(col="outcome_thrombo_50",   label="Thrombopénie < 50 G/L"),
  list(col="outcome_composite",    label="Outcome composite")
)

ate_secondary <- rbindlist(lapply(sec_outcomes, function(oc) {
  run_aipw_with_ci(oc$col, oc$label)
}))

# Sauvegarder
safe_write_csv(ate_secondary, file.path(TABLES_DIR, "ate_results_secondary_dose_based.csv"), overwrite=TRUE)
log_msg("Outcomes secondaires dose-based sauvegardés")

# Tableau ATE complet (dose-based principale + secondaires + GRF)
ate_complete <- rbind(ate_updated, ate_secondary, fill=TRUE)
safe_write_csv(ate_complete, file.path(TABLES_DIR, "ate_results_complete.csv"), overwrite=TRUE)
log_msg("Tableau ATE complet sauvegardé")
print(ate_complete[, .(estimator, outcome, estimate_rd, ci_lower, ci_upper)])

timer_end(t0, "A4+A5+A11 outcomes")
log_msg("=== STARTING ===")
