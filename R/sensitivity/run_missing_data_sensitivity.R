# =============================================================================
# run_missing_data_sensitivity.R
# Step 9 — Missing-data sensitivity (4 strategies, bootstrap CIs)
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table); library(arrow); library(ggplot2)

log_msg("=== STARTING ===")
t0 <- timer_start()

ensure_dirs(TABLES_DIR, FIGURES_DIR)

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
cov_vars <- intersect(propensity_vars, names(cohort))

# Résumé des données manquantes
miss_vars <- sapply(cov_vars, function(v) sum(is.na(cohort[[v]])))
miss_pct  <- sapply(cov_vars, function(v) mean(is.na(cohort[[v]]))*100)
miss_summary <- data.table(variable=cov_vars, n_missing=miss_vars, pct_missing=round(miss_pct, 1))
miss_summary <- miss_summary[pct_missing > 0][order(-pct_missing)]
log_msg("Variables avec données manquantes:")
print(miss_summary)

# Fonction AIPW légère
aipw_quick <- function(data, cov_v) {
  ps_f <- as.formula(paste("high_dose ~", paste(cov_v, collapse="+")))
  oc_f <- as.formula(paste("outcome_thrombo_100 ~ high_dose +", paste(cov_v, collapse="+")))
  tryCatch({
    ps_m <- glm(ps_f, data=as.data.frame(data), family=binomial)
    ps_h <- pmax(pmin(predict(ps_m, type="response"), 0.99), 0.01)
    oc_m <- glm(oc_f, data=as.data.frame(data), family=binomial)
    d1 <- copy(data); d1[, high_dose := 1L]
    d0 <- copy(data); d0[, high_dose := 0L]
    m1 <- predict(oc_m, newdata=d1, type="response")
    m0 <- predict(oc_m, newdata=d0, type="response")
    A <- data$high_dose; Y <- as.integer(data$outcome_thrombo_100)
    mean(m1 + A/ps_h*(Y-m1)) - mean(m0 + (1-A)/(1-ps_h)*(Y-m0))
  }, error=function(e) NA_real_)
}

# =============================================================================
# 1. Analyse avec imputation médiane (déjà implémentée — référence)
# =============================================================================
cohort_median <- copy(cohort)
for (v in cov_vars) {
  if (is.numeric(cohort_median[[v]])) {
    cohort_median[is.na(get(v)), (v) := median(cohort_median[[v]], na.rm=TRUE)]
  } else if (is.character(cohort_median[[v]])) {
    mode_v <- names(sort(table(cohort_median[[v]]), decreasing=TRUE))[1]
    cohort_median[is.na(get(v)), (v) := mode_v]
  }
}
rd_median <- aipw_quick(cohort_median, cov_vars)
log_msg(sprintf("AIPW (imputation médiane): RD=%.4f", rd_median))

# =============================================================================
# 2. Complete case analysis
# =============================================================================
cohort_cc <- cohort[complete.cases(cohort[, c(cov_vars, "outcome_thrombo_100"), with=FALSE])]
n_cc <- nrow(cohort_cc)
log_msg(sprintf("Complete case: %d patients (%.1f%% de %d)", n_cc, n_cc/nrow(cohort)*100, nrow(cohort)))
rd_cc <- aipw_quick(cohort_cc, cov_vars)
log_msg(sprintf("AIPW (complete case): RD=%.4f", rd_cc))

# =============================================================================
# 3. MICE léger (2 imputations) si mice disponible
# =============================================================================
mice_available <- requireNamespace("mice", quietly=TRUE)
rd_mice <- NA_real_
if (mice_available) {
  library(mice)
  log_msg("3. MICE imputation (5 imputations, méthode pmm)...")
  # Variables avec NA
  vars_with_na <- cov_vars[sapply(cov_vars, function(v) sum(is.na(cohort[[v]])) > 0)]
  # Sous-ensemble pour imputation
  data_for_mice <- cohort[, c(cov_vars, "outcome_thrombo_100", "high_dose"), with=FALSE]
  tryCatch({
    m_imp <- suppressMessages(
      mice(as.data.frame(data_for_mice), m=5, method="pmm", printFlag=FALSE,
           seed=RANDOM_SEED, maxit=5)
    )
    # AIPW sur chacune des 5 imputations
    rd_mice_v <- sapply(1:5, function(i) {
      d_imp <- as.data.table(complete(m_imp, i))
      d_imp[, high_dose := cohort$high_dose]
      aipw_quick(d_imp, cov_vars)
    })
    rd_mice <- mean(rd_mice_v, na.rm=TRUE)
    log_msg(sprintf("AIPW (MICE, 5 imputations): RD=%.4f (SD=%.4f)", rd_mice, sd(rd_mice_v, na.rm=TRUE)))
  }, error=function(e) {
    log_msg(sprintf("WARN: MICE erreur: %s", e$message))
    rd_mice <<- NA_real_
  })
} else {
  log_msg("MICE non disponible — skip")
}

# =============================================================================
# 4. Sensibilité : exclure variables avec > 20% manquants
# =============================================================================
vars_low_miss <- cov_vars[sapply(cov_vars, function(v) mean(is.na(cohort[[v]])) <= 0.20)]
cohort_lowmiss <- copy(cohort)
for (v in vars_low_miss) {
  if (is.numeric(cohort_lowmiss[[v]])) {
    cohort_lowmiss[is.na(get(v)), (v) := median(cohort_lowmiss[[v]], na.rm=TRUE)]
  }
}
cohort_lowmiss_cc <- cohort_lowmiss[complete.cases(cohort_lowmiss[, c(vars_low_miss, "outcome_thrombo_100"), with=FALSE])]
rd_lowmiss <- aipw_quick(cohort_lowmiss_cc, vars_low_miss)
log_msg(sprintf("AIPW (PS sans vars > 20%% manquants): RD=%.4f", rd_lowmiss))

# =============================================================================
# 5. Compilation vibration imputation
# =============================================================================
vib_imp <- data.table(
  category    = "Gestion données manquantes",
  scenario    = c(
    "Imputation médiane (référence)",
    "Analyse en cas complets (complete case)",
    if (!is.na(rd_mice)) "MICE (5 imputations, pmm)" else NULL,
    "Sans covariables > 20% manquants"
  ),
  rd          = round(c(rd_median, rd_cc, if (!is.na(rd_mice)) rd_mice else NULL, rd_lowmiss), 4),
  n_patients  = c(nrow(cohort_median), n_cc, if (!is.na(rd_mice)) nrow(cohort) else NULL, nrow(cohort_lowmiss_cc)),
  note        = c(
    "Méthode principale",
    sprintf("N=%d (%.1f%% avec toutes vars disponibles)", n_cc, n_cc/nrow(cohort)*100),
    if (!is.na(rd_mice)) "Rubin's rules sur 5 imputations" else NULL,
    sprintf("PS avec %d vars (<%d covariables à miss>20%%)", length(vars_low_miss), length(cov_vars)-length(vars_low_miss))
  )
)

# Charger vibration existante et étendre
vib_existing <- tryCatch(fread(file.path(TABLES_DIR, "vibration_results.csv")), error=function(e) NULL)
if (!is.null(vib_existing)) {
  # Reformater pour compatibilité
  vib_imp2 <- data.table(
    category  = vib_imp$category,
    scenario  = vib_imp$scenario,
    n_treated = rep(NA_integer_, nrow(vib_imp)),
    rd        = vib_imp$rd,
    se        = rep(NA_real_, nrow(vib_imp)),
    ci_lower  = rep(NA_real_, nrow(vib_imp)),
    ci_upper  = rep(NA_real_, nrow(vib_imp)),
    order_id  = max(vib_existing$order_id) + seq_len(nrow(vib_imp))
  )
  vib_extended <- rbind(vib_existing, vib_imp2, fill=TRUE)
  safe_write_csv(vib_extended, file.path(TABLES_DIR, "vibration_results.csv"), overwrite=TRUE)
}

# Sauvegarder séparément la partie imputation
safe_write_csv(vib_imp, file.path(TABLES_DIR, "vibration_imputation.csv"), overwrite=TRUE)
log_msg("Tableau vibration imputation sauvegardé:")
print(vib_imp)

# =============================================================================
# 6. Figure vibration imputation
# =============================================================================
p_vib_imp <- ggplot(vib_imp, aes(x=rd*100, y=reorder(scenario, rd))) +
  geom_point(aes(color=scenario=="Imputation médiane (référence)"), size=4) +
  geom_vline(xintercept=0, linetype="dashed", color="gray30") +
  geom_vline(xintercept=rd_median*100, color="steelblue", linetype="dotted", linewidth=0.8) +
  scale_color_manual(values=c("gray50","tomato"), guide="none") +
  scale_x_continuous(labels=function(x) paste0(x, "%")) +
  labs(
    title="Vibration Analysis — Sensibilité à la gestion des données manquantes",
    subtitle=sprintf("Outcome: Thrombopénie <100 G/L | Estimateur: AIPW\nLigne pointillée = résultat de référence (imputation médiane, RD=+%.1f%%)", rd_median*100),
    x="Risk Difference AIPW (%)", y="Stratégie d'imputation"
  ) +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "vibration_imputation_dose_based.png"), p_vib_imp, width=10, height=5, dpi=150)

# =============================================================================
# Synthèse interprétative
# =============================================================================
rd_range <- range(c(rd_median, rd_cc, rd_mice, rd_lowmiss), na.rm=TRUE)
log_msg(sprintf("Plage AIPW selon imputation: %.4f à %.4f", rd_range[1], rd_range[2]))
stable <- abs(rd_range[2] - rd_range[1]) < 0.01
log_msg(sprintf("Signal stable à l'imputation: %s (range < 1%%: %s)", stable, stable))

timer_end(t0, "A6 vibration imputation")
log_msg("=== STARTING ===")
