# =============================================================================
# run_grf_calibration.R
# Step 10 — GRF calibration tests and variable importance
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table); library(arrow); library(grf); library(ggplot2)

log_msg("=== STARTING ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR, DIAGNOSTICS_DIR, MODELS_DIR)

# Charger modèle CF existant
cf <- readRDS(file.path(MODELS_DIR, "causal_forest.rds"))
log_msg(sprintf("Causal Forest chargé: %d arbres", cf$`_num_trees`))

# =============================================================================
# 1. test_calibration() — test global d'hétérogénéité
# =============================================================================
log_msg("1. test_calibration()...")
calib_test <- test_calibration(cf)
log_msg("Résultat test_calibration:")
print(calib_test)

# test_calibration retourne un objet 'coeftest' (matrice numérique avec dimnames)
# Structure: 2 lignes x 4 colonnes [Estimate, Std.Error, t value, Pr(>t)]
calib_mat <- as.matrix(calib_test)  # déjà une matrice numérique
mean_est  <- calib_mat["mean.forest.prediction", "Estimate"]
mean_p    <- calib_mat["mean.forest.prediction", "Pr(>t)"]
diff_est  <- calib_mat["differential.forest.prediction", "Estimate"]
diff_p    <- calib_mat["differential.forest.prediction", "Pr(>t)"]

log_msg(sprintf("Calibration: mean_est=%.3f, mean_p=%.4f, diff_est=%.3f, diff_p=%.4f",
                mean_est, mean_p, diff_est, diff_p))

log_msg(sprintf("Mean forest prediction: est=%.3f, p=%.4f", mean_est, mean_p))
log_msg(sprintf("Differential forest prediction: est=%.3f, p=%.4f", diff_est, diff_p))

# =============================================================================
# 2. ATE via GRF (pour intégration dans le tableau ATE principal)
# =============================================================================
log_msg("2. ATE via GRF...")
ate_grf <- average_treatment_effect(cf, target.sample="all")
ate_grf_dt <- data.table(
  test     = "ATE (GRF, target.sample=all)",
  estimate = round(ate_grf["estimate"], 4),
  std_err  = round(ate_grf["std.err"], 4),
  ci_lower = round(ate_grf["estimate"] - 1.96*ate_grf["std.err"], 4),
  ci_upper = round(ate_grf["estimate"] + 1.96*ate_grf["std.err"], 4)
)
log_msg(sprintf("ATE GRF: %.4f [%.4f, %.4f]",
                ate_grf["estimate"],
                ate_grf["estimate"] - 1.96*ate_grf["std.err"],
                ate_grf["estimate"] + 1.96*ate_grf["std.err"]))

# =============================================================================
# 3. Calibration : observed vs predicted CATE
# =============================================================================
log_msg("3. Plot calibration observed vs predicted...")

cate_pred <- predict(cf)$predictions
cohort_clean <- as.data.table(read_parquet(file.path(DATA_FINAL, "features_cohort.parquet")))
cohort_clean[, cate_hat := cate_pred]
cohort_clean[, cate_decile := cut(cate_hat,
  breaks = quantile(cate_hat, seq(0, 1, 0.1), na.rm=TRUE),
  include.lowest = TRUE, labels=FALSE
)]

# Pour chaque décile de CATE : comparer ATE observé (AIPW simple)
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
cov_vars <- intersect(propensity_vars, names(cohort_clean))
for (v in cov_vars) {
  if (is.numeric(cohort_clean[[v]])) cohort_clean[is.na(get(v)), (v) := median(cohort_clean[[v]], na.rm=TRUE)]
}

ps_formula <- as.formula(paste("high_dose ~", paste(cov_vars, collapse="+")))
ps_m <- glm(ps_formula, data=as.data.frame(cohort_clean), family=binomial)
ps_hat <- predict(ps_m, type="response")

decile_stats <- cohort_clean[!is.na(cate_decile), .(
  cate_mean  = mean(cate_hat),
  n          = .N,
  crude_rd   = mean(outcome_thrombo_100[high_dose==1]) - mean(outcome_thrombo_100[high_dose==0]),
  n_treated  = sum(high_dose),
  n_control  = sum(1-high_dose)
), by=cate_decile][order(cate_decile)]

# Calibration figure
p_calib <- ggplot(decile_stats[n_treated >= 5 & n_control >= 5],
                   aes(x=cate_mean*100, y=crude_rd*100)) +
  geom_point(aes(size=n), color="steelblue", alpha=0.8) +
  geom_smooth(method="lm", se=TRUE, color="darkorange", linewidth=0.8) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
  geom_hline(yintercept=0, color="gray70", linetype="dotted") +
  geom_vline(xintercept=0, color="gray70", linetype="dotted") +
  scale_size_continuous(range=c(2,8), name="N patients") +
  labs(
    title="Calibration du Causal Forest (GRF)",
    subtitle=sprintf("CATE prédit (axe X) vs RD observé par décile (axe Y)\nTest calibration : mean_pred p=%.3f | diff_pred p=%.3f",
                     mean_p, diff_p),
    x="CATE moyen prédit (%)", y="RD observé (brut, par décile) (%)"
  ) +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "grf_calibration_plot.png"), p_calib, width=9, height=6, dpi=150)

# =============================================================================
# 4. Décision sur l'hétérogénéité et les variable importance plots
# =============================================================================
hte_evidence <- if (!is.na(diff_p) && diff_p < 0.05) {
  "DÉTECTÉE (différential forest prediction p < 0.05)"
} else if (!is.na(diff_p) && diff_p < 0.1) {
  "MARGINALE (différential forest prediction 0.05 ≤ p < 0.10)"
} else {
  "NON DÉMONTRÉE (différential forest prediction p ≥ 0.10)"
}
log_msg(sprintf("Hétérogénéité GRF: %s", hte_evidence))

# =============================================================================
# 5. Tableau synthèse calibration
# =============================================================================
calib_summary <- data.table(
  test = c("mean.forest.prediction", "differential.forest.prediction"),
  description = c(
    "Test si ATE moyen prédit ≈ ATE observé (calibration globale)",
    "Test si hétérogénéité du CATE est supérieure à un modèle constant (test HTE)"
  ),
  estimate = round(c(mean_est, diff_est), 4),
  p_value  = round(c(mean_p, diff_p), 4),
  interpretation = c(
    ifelse(mean_p < 0.05, "Calibration globale adéquate", "Mauvaise calibration globale"),
    hte_evidence
  )
)
safe_write_csv(calib_summary, file.path(TABLES_DIR, "grf_calibration_tests.csv"), overwrite=TRUE)
log_msg("Tableau calibration sauvegardé")

# ATE GRF pour intégration dans tableau principal
safe_write_csv(ate_grf_dt, file.path(TABLES_DIR, "grf_ate.csv"), overwrite=TRUE)

# =============================================================================
# 6. Figure variable importance — conservée AVEC mention explicite des limites
# =============================================================================
var_imp_dt <- fread(file.path(TABLES_DIR, "grf_variable_importance.csv"))
hte_caveat <- if (diff_p >= 0.05) {
  "ATTENTION : hétérogénéité non démontrée formellement — importance à interpréter avec prudence"
} else {
  "Hétérogénéité démontrée formellement (p < 0.05)"
}

p_varimp <- ggplot(head(var_imp_dt[!is.na(importance) & importance > 0], 15),
                    aes(x=importance, y=reorder(variable, importance))) +
  geom_col(fill="steelblue", alpha=0.75) +
  labs(
    title="Variables associées à l'hétérogénéité (Causal Forest GRF)",
    subtitle=sprintf("%s\np(diff.forest)=%.3f | ATE GRF=%.3f [%.3f, %.3f]",
                     hte_caveat, diff_p,
                     ate_grf["estimate"],
                     ate_grf["estimate"] - 1.96*ate_grf["std.err"],
                     ate_grf["estimate"] + 1.96*ate_grf["std.err"]),
    x="Importance relative (GRF)", y="Variable"
  ) +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "hte_variable_importance.png"), p_varimp, width=10, height=7, dpi=150)

# =============================================================================
# Rapport
# =============================================================================
report_a7 <- sprintf(
'# Annexe A7 — Test de Calibration et Hétérogénéité (GRF)
**Date :** %s

## Méthode

La fonction `test_calibration()` du package grf (Tibshirani et al.) évalue :
1. **mean.forest.prediction** : si la prédiction CATE moyenne correspond à l\'ATE observé (calibration globale)
2. **differential.forest.prediction** : si le CATE varie au-delà d\'un modèle constant (test d\'hétérogénéité)

Un test différentiel significatif est requis pour justifier l\'interprétation des importances de variables.

## Résultats

| Test | Estimate | p-value | Interprétation |
|------|----------|---------|----------------|
| mean.forest.prediction | %.3f | %.4f | %s |
| differential.forest.prediction | %.3f | %.4f | %s |

## Conclusion

**Hétérogénéité : %s**

%s

## ATE global via GRF

ATE = %.4f ± %.4f [IC 95%% : %.4f, %.4f]

## Implications pour le rapport

%s
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  mean_est, mean_p, calib_summary$interpretation[1],
  diff_est, diff_p, hte_evidence,
  hte_evidence,
  if (diff_p >= 0.05)
    "Le test différentiel n'est pas significatif. Les variable importance plots doivent être interprétés avec prudence — ils indiquent quelles variables *pourraient* modérer l'effet si une HTE existait, mais leur pertinence formelle n'est pas établie dans ces données."
  else
    "Le test différentiel est significatif. Les variable importance plots peuvent être interprétés, avec prudence.",
  ate_grf["estimate"], ate_grf["std.err"],
  ate_grf["estimate"] - 1.96*ate_grf["std.err"],
  ate_grf["estimate"] + 1.96*ate_grf["std.err"],
  if (diff_p >= 0.05)
    paste("- Les phrases affirmant que telle variable 'drive l'hétérogénéité' sont retirées du message principal\n",
          "- Les figures d'importance sont déplacées en exploratoire/annexe avec caveat explicite\n",
          "- L'ATE GRF reste un résultat valide et est intégré au tableau principal")
  else
    "- Les résultats HTE peuvent être présentés avec une interprétation prudente"
)
writeLines(report_a7, file.path(REPORTS_DIR, "appendix_A7_grf_calibration.md"))
log_msg(sprintf("Rapport A7 écrit | diff_p=%.4f | %s", diff_p, hte_evidence))

timer_end(t0, "A7 GRF calibration")
log_msg("=== STARTING ===")
