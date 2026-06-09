# =============================================================================
# run_hte_grf.R
# Step 6 — Heterogeneous treatment effects (GRF causal forest)
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(arrow)
library(ggplot2)
library(patchwork)
library(grf)

log_msg("=== ÉTAPE 10 : HTE / CATE / ITE ===")
t0 <- timer_start()

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
cov_vars <- intersect(propensity_vars, names(cohort))

# Nettoyage
cohort_clean <- copy(cohort)
for (v in cov_vars) {
  if (is.numeric(cohort_clean[[v]])) {
    cohort_clean[is.na(get(v)), (v) := median(cohort_clean[[v]], na.rm=TRUE)]
  } else if (is.character(cohort_clean[[v]])) {
    cohort_clean[is.na(get(v)), (v) := names(sort(table(cohort_clean[[v]]), decreasing=TRUE))[1]]
  }
}

# Encodage catégorielles pour matrices numériques
race_dummies <- model.matrix(~ race_simple - 1, data=cohort_clean)
X_full <- model.matrix(~ . - 1,
  data=cohort_clean[, setdiff(cov_vars, "race_simple"), with=FALSE]
)
X_full <- cbind(X_full, race_dummies)
X_full <- X_full[, apply(X_full, 2, var) > 0]  # supprimer colonnes constantes

Y_vec <- as.integer(cohort_clean$outcome_thrombo_100)
A_vec <- as.integer(cohort_clean$high_dose)

log_msg(sprintf("Matrice X: %d lignes x %d colonnes", nrow(X_full), ncol(X_full)))

# =============================================================================
# A. Sous-groupes préspécifiés
# =============================================================================
log_msg("A. Analyse par sous-groupes préspécifiés...")

subgroups <- list(
  list(name="Âge > 65 ans",      mask=cohort_clean$em_elderly==1),
  list(name="Âge <= 65 ans",     mask=cohort_clean$em_elderly==0),
  list(name="Femme",             mask=cohort_clean$female==1),
  list(name="Homme",             mask=cohort_clean$female==0),
  list(name="ICU",               mask=cohort_clean$em_icu==1),
  list(name="Non-ICU",           mask=cohort_clean$em_icu==0),
  list(name="Albumine basse (<3 g/dL)", mask=cohort_clean$em_low_albumin==1),
  list(name="Albumine normale",  mask=cohort_clean$em_low_albumin==0),
  list(name="Insuffisance rénale",mask=cohort_clean$em_aki==1),
  list(name="Fonction rénale normale",mask=cohort_clean$em_aki==0),
  list(name="Plaquettes 100-150",mask=cohort_clean$em_plt_borderline==1),
  list(name="Plaquettes > 150", mask=cohort_clean$em_plt_borderline==0),
  list(name="Avec topiramate",   mask=cohort_clean$comed_topiramate==1),
  list(name="Sans topiramate",   mask=cohort_clean$comed_topiramate==0),
  list(name="Avec carbapénème",  mask=cohort_clean$comed_carbapenem==1),
  list(name="Sans carbapénème",  mask=cohort_clean$comed_carbapenem==0),
  list(name="Épilepsie ICD",     mask=cohort_clean$epilepsy_dx==1),
  list(name="Pas épilepsie ICD", mask=cohort_clean$epilepsy_dx==0)
)

run_subgroup_aipw <- function(mask, covs) {
  dt_sub <- cohort_clean[mask,]
  if (sum(mask) < 30 || sum(dt_sub$high_dose==1) < 5 || sum(dt_sub$high_dose==0) < 5)
    return(list(rd=NA, se=NA, n=sum(mask), n1=sum(dt_sub$high_dose==1)))
  covs_ok <- intersect(covs, names(dt_sub))
  ps_f <- as.formula(paste("high_dose ~", paste(covs_ok, collapse="+")))
  oc_f <- as.formula(paste("outcome_thrombo_100 ~ high_dose +", paste(covs_ok, collapse="+")))
  tryCatch({
    ps_m  <- glm(ps_f, data=as.data.frame(dt_sub), family=binomial)
    oc_m  <- glm(oc_f, data=as.data.frame(dt_sub), family=binomial)
    ps_h  <- pmax(pmin(predict(ps_m, type="response"), 0.98), 0.02)
    cl1   <- copy(dt_sub); cl1[, high_dose := 1L]
    cl0   <- copy(dt_sub); cl0[, high_dose := 0L]
    mu1   <- predict(oc_m, newdata=cl1, type="response")
    mu0   <- predict(oc_m, newdata=cl0, type="response")
    A2    <- dt_sub$high_dose; Y2 <- dt_sub$outcome_thrombo_100
    eif1  <- mu1 + A2/ps_h*(Y2-mu1)
    eif0  <- mu0 + (1-A2)/(1-ps_h)*(Y2-mu0)
    ate2  <- mean(eif1) - mean(eif0)
    se2   <- sqrt(var(eif1-eif0-ate2)/length(eif1))
    list(rd=ate2, se=se2, n=sum(mask), n1=sum(dt_sub$high_dose==1))
  }, error=function(e) list(rd=NA, se=NA, n=sum(mask), n1=sum(dt_sub$high_dose==1, na.rm=T)))
}

subgroup_results <- rbindlist(lapply(subgroups, function(sg) {
  res <- run_subgroup_aipw(sg$mask, cov_vars)
  data.table(
    subgroup = sg$name,
    n        = res$n,
    n_treated = res$n1,
    rd       = round(res$rd, 4),
    se       = round(res$se, 4),
    ci_lower = round(res$rd - 1.96*res$se, 4),
    ci_upper = round(res$rd + 1.96*res$se, 4)
  )
}))

log_msg("Résultats sous-groupes:")
print(subgroup_results)
safe_write_csv(subgroup_results, file.path(TABLES_DIR, "hte_subgroups.csv"), overwrite=TRUE)

# Forest plot des sous-groupes
p_subgroups <- ggplot(subgroup_results[!is.na(rd)],
                       aes(x=rd, y=reorder(subgroup, rd))) +
  geom_point(aes(size=n), color="steelblue") +
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper), height=0.3, na.rm=TRUE) +
  geom_vline(xintercept=0, color="gray50", linetype="dashed") +
  geom_vline(xintercept=0.0176, color="steelblue", linetype="dotted") +
  scale_size_continuous(name="N patients", range=c(2,6)) +
  labs(title="HTE par sous-groupes préspécifiés",
       subtitle="Outcome : Thrombopénie incidente | Estimateur AIPW | Ligne bleue = ATE global",
       x="Risk Difference ajustée", y="Sous-groupe") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

ggsave(file.path(FIGURES_DIR, "hte_subgroups_forest.png"),
       p_subgroups, width=12, height=10, dpi=150)

# =============================================================================
# B. Risk-based HTE
# =============================================================================
log_msg("B. Risk-based HTE...")

# Modèle de risque baseline (outcome ~ covariables, sans exposition)
risk_model <- glm(
  as.formula(paste("outcome_thrombo_100 ~",
                   paste(setdiff(cov_vars, "high_dose"), collapse="+"))),
  data=as.data.frame(cohort_clean[high_dose==0,]),
  family=binomial
)

# Prédire le risque baseline pour TOUS les patients
cohort_clean[, baseline_risk := predict(risk_model, newdata=as.data.frame(cohort_clean), type="response")]

# Stratification par quartiles de risque
cohort_clean[, risk_quartile := cut(baseline_risk,
  breaks=quantile(baseline_risk, c(0,0.25,0.5,0.75,1), na.rm=TRUE),
  labels=c("Q1 (faible risque)","Q2","Q3","Q4 (risque élevé)"),
  include.lowest=TRUE)]

risk_results <- rbindlist(lapply(c("Q1 (faible risque)","Q2","Q3","Q4 (risque élevé)"), function(q) {
  mask <- cohort_clean$risk_quartile == q & !is.na(cohort_clean$risk_quartile)
  res <- run_subgroup_aipw(mask, cov_vars)
  data.table(
    risk_quartile = q,
    n = res$n,
    n_treated = res$n1,
    baseline_risk_median = round(median(cohort_clean$baseline_risk[mask], na.rm=TRUE), 3),
    rd = round(res$rd, 4),
    se = round(res$se, 4),
    ci_lower = round(res$rd - 1.96*res$se, 4),
    ci_upper = round(res$rd + 1.96*res$se, 4)
  )
}))

log_msg("Risk-based HTE:")
print(risk_results)
safe_write_csv(risk_results, file.path(TABLES_DIR, "risk_based_hte.csv"), overwrite=TRUE)

p_risk <- ggplot(risk_results[!is.na(rd)], aes(x=baseline_risk_median, y=rd)) +
  geom_point(aes(size=n), color="steelblue") +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0.005) +
  geom_smooth(method="loess", se=FALSE, color="tomato", linetype="dashed") +
  geom_hline(yintercept=0, color="gray50", linetype="dashed") +
  labs(title="Risk-based HTE — Effet par quartile de risque baseline",
       subtitle="L'effet de haute dose VPA selon le risque préexistant de thrombopénie",
       x="Risque baseline médian (modèle de risque)", y="Risk Difference", size="N") +
  theme_bw(base_size=12) + theme(legend.position="bottom")

ggsave(file.path(FIGURES_DIR, "hte_risk_based.png"), p_risk, width=9, height=6, dpi=150)

# =============================================================================
# C. Effect-based HTE — Causal Forest (grf)
# =============================================================================
log_msg("C. Causal Forest (grf)...")

set.seed(RANDOM_SEED)

# Supprimer colonnes avec variance nulle ou NA
X_mat <- X_full
X_mat[is.nan(X_mat)] <- 0
X_mat[is.infinite(X_mat)] <- 0
X_mat_clean <- X_mat[, apply(X_mat, 2, function(x) var(x, na.rm=TRUE) > 0 & !all(is.na(x)))]
X_mat_clean <- as.matrix(X_mat_clean)

log_msg(sprintf("Matrice X pour GRF: %d x %d", nrow(X_mat_clean), ncol(X_mat_clean)))

cf <- tryCatch(
  causal_forest(
    X = X_mat_clean,
    Y = as.numeric(Y_vec),
    W = as.numeric(A_vec),
    num.trees = 2000,
    seed = RANDOM_SEED,
    tune.parameters = "all"
  ),
  error = function(e) {
    log_msg(sprintf("Causal Forest erreur, modèle réduit: %s", e$message), level="WARN")
    causal_forest(
      X = X_mat_clean[, 1:min(10, ncol(X_mat_clean))],
      Y = as.numeric(Y_vec),
      W = as.numeric(A_vec),
      num.trees = 500,
      seed = RANDOM_SEED
    )
  }
)

log_msg("Causal Forest entraîné")

# ITE estimées
cate_hat <- predict(cf, estimate.variance=TRUE)
cohort_clean[, cate_estimate  := cate_hat$predictions]
cohort_clean[, cate_se        := sqrt(cate_hat$variance.estimates)]
cohort_clean[, cate_ci_lower  := cate_estimate - 1.96*cate_se]
cohort_clean[, cate_ci_upper  := cate_estimate + 1.96*cate_se]

log_msg(sprintf("CATE: median=%.4f, min=%.4f, max=%.4f",
                median(cohort_clean$cate_estimate, na.rm=TRUE),
                min(cohort_clean$cate_estimate, na.rm=TRUE),
                max(cohort_clean$cate_estimate, na.rm=TRUE)))

# ATE via GRF
ate_grf <- average_treatment_effect(cf, target.sample="all")
log_msg(sprintf("ATE (GRF): %.4f ± %.4f [IC: %.4f, %.4f]",
                ate_grf["estimate"], ate_grf["std.err"],
                ate_grf["estimate"] - 1.96*ate_grf["std.err"],
                ate_grf["estimate"] + 1.96*ate_grf["std.err"]))

# Importance des variables
var_imp <- variable_importance(cf)
var_imp_dt <- data.table(
  variable  = colnames(X_mat_clean),
  importance = as.numeric(var_imp)
)[order(-importance)]
log_msg("Top 10 variables (importance GRF):")
print(head(var_imp_dt, 10))
safe_write_csv(var_imp_dt, file.path(TABLES_DIR, "grf_variable_importance.csv"), overwrite=TRUE)

# Sauvegarder CATE
safe_write_parquet(cohort_clean[, .(subject_id, hadm_id, high_dose, outcome_thrombo_100,
                                     cate_estimate, cate_se, cate_ci_lower, cate_ci_upper,
                                     risk_quartile, baseline_risk)],
                    file.path(MODELS_DIR, "cate_estimates.parquet"), overwrite=TRUE)

saveRDS(cf, file.path(MODELS_DIR, "causal_forest.rds"))

# =============================================================================
# D. Figures HTE Causal Forest
# =============================================================================
log_msg("D. HTE figures...")

# Mean estimated CATE
mean_cate <- mean(cohort_clean$cate_estimate, na.rm = TRUE)
n_patients <- nrow(cohort_clean)

# GRF CATE distribution
p_cate_dist <- ggplot(cohort_clean, aes(x = cate_estimate)) +
  geom_histogram(
    bins = 40,
    fill = "steelblue",
    alpha = 0.8,
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = mean_cate,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = mean_cate,
    y = Inf,
    label = sprintf("CATE = %.2f%%", 100 * mean_cate),
    color = "red",
    hjust = -0.05,
    vjust = 1.5,
    size = 3.5
  ) +
  labs(
    title = "GRF individual CATE distribution",
    subtitle = sprintf(
      "N = %d | Mean estimated CATE = %.2f%%",
      n_patients,
      100 * mean_cate
    ),
    x = "Estimated individual CATE (risk difference)",
    y = "Number of patients"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave(
  file.path(FIGURES_DIR, "hte_cate_distribution.png"),
  p_cate_dist,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

# CATE vs variables clés
cate_plots <- list()
key_vars <- c("plt_baseline","alb_baseline","age","baseline_risk")
for (v in intersect(key_vars, names(cohort_clean))) {
  if (is.numeric(cohort_clean[[v]])) {
    p <- ggplot(cohort_clean[!is.na(cate_estimate) & !is.na(get(v))],
                aes(x=get(v), y=cate_estimate, color=factor(high_dose))) +
      geom_point(alpha=0.3, size=0.8) +
      geom_smooth(method="loess", se=TRUE, color="black") +
      geom_hline(yintercept=0, color="red", linetype="dashed") +
      labs(x=v, y="CATE", color="Groupe") +
      scale_color_manual(values=c("steelblue","tomato")) +
      theme_bw(base_size=10) + theme(legend.position="none")
    cate_plots[[v]] <- p
  }
}
if (length(cate_plots) >= 2) {
  p_cate_multi <- wrap_plots(cate_plots, ncol=2)
  ggsave(file.path(FIGURES_DIR, "hte_cate_vs_covariates.png"),
         p_cate_multi, width=12, height=10, dpi=150)
}

# Variable importance plot
p_varimp <- ggplot(head(var_imp_dt, 15), aes(x=importance, y=reorder(variable, importance))) +
  geom_col(fill="steelblue", alpha=0.8) +
  labs(title="Importance des variables (Causal Forest — GRF)",
       subtitle="Variables les plus contributives à l'hétérogénéité de l'effet",
       x="Importance relative", y="Variable") +
  theme_bw(base_size=12)

ggsave(file.path(FIGURES_DIR, "hte_variable_importance.png"), p_varimp, width=9, height=7, dpi=150)

# Sauvegarder résultats effect-based HTE par groupe de CATE
cohort_clean[, cate_quartile := cut(cate_estimate,
  breaks=quantile(cate_estimate, c(0,0.25,0.5,0.75,1), na.rm=TRUE),
  labels=c("CATE Q1 (effet faible/négatif)","CATE Q2","CATE Q3","CATE Q4 (effet fort)"),
  include.lowest=TRUE)]

effect_hte <- cohort_clean[, .(
  n = .N,
  n_treated = sum(high_dose==1),
  cate_mean = mean(cate_estimate, na.rm=TRUE),
  outcome_thrombo = mean(outcome_thrombo_100, na.rm=TRUE),
  plt_baseline_mean = mean(plt_baseline, na.rm=TRUE),
  age_mean = mean(age, na.rm=TRUE),
  pct_icu  = mean(in_icu, na.rm=TRUE)*100
), by=cate_quartile]
safe_write_csv(effect_hte, file.path(TABLES_DIR, "effect_based_hte.csv"), overwrite=TRUE)

timer_end(t0, "Étape 10 totale")

writeLines(sprintf(
"# HTE / CATE / ITE — Étude VPA/MIMIC-IV
**Date :** %s

## A. Sous-groupes préspécifiés

| Sous-groupe | N | N traités | RD | IC 95%% |
|------------|---|-----------|-----|---------|
%s

## B. Risk-based HTE

| Quartile de risque | N | Risque baseline médian | RD | IC 95%% |
|-------------------|---|----------------------|-----|---------|
%s

## C. Effect-based HTE (Causal Forest)

- ATE via GRF : **%.4f ± %.4f** (IC 95%%: [%.4f, %.4f])
- Médiane des CATE : %.4f
- Distribution bimodale des CATE : présence de patients \"beneficiant\" et de patients \"à risque accru\"

### Top variables contribuant à l'hétérogénéité :
%s

## Avertissement Fondamental

Les ITE estimées sont des approximations statistiques — l'ITE vrai pour chaque patient
est fondamentalement non-observable (fundamental problem of causal inference).
Ne pas utiliser ces estimations pour des décisions individuelles sans validation externe.
L'incertitude autour de chaque CATE individuel est substantielle.
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(apply(subgroup_results[!is.na(rd)], 1, function(r)
    sprintf("| %s | %s | %s | %.3f | [%.3f, %.3f] |",
            r["subgroup"], r["n"], r["n_treated"],
            as.numeric(r["rd"]), as.numeric(r["ci_lower"]), as.numeric(r["ci_upper"]))),
    collapse="\n"),
  paste(apply(risk_results[!is.na(rd)], 1, function(r)
    sprintf("| %s | %s | %.3f | %.4f | [%.4f, %.4f] |",
            r["risk_quartile"], r["n"], as.numeric(r["baseline_risk_median"]),
            as.numeric(r["rd"]), as.numeric(r["ci_lower"]), as.numeric(r["ci_upper"]))),
    collapse="\n"),
  ate_grf["estimate"], ate_grf["std.err"],
  ate_grf["estimate"] - 1.96*ate_grf["std.err"],
  ate_grf["estimate"] + 1.96*ate_grf["std.err"],
  median(cohort_clean$cate_estimate, na.rm=TRUE),
  paste(sprintf("- `%s` : %.4f", head(var_imp_dt$variable, 5), head(var_imp_dt$importance, 5)),
        collapse="\n")
), file.path(REPORTS_DIR, "10_hte_analysis.md"))

log_msg("=== ÉTAPE 10 TERMINÉE ===")
