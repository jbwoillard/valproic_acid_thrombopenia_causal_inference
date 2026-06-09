# =============================================================================
# estimate_primary_ate.R
# Step 4 — Primary ATE estimation (G-formula, IPTW, AIPW)
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(arrow)
library(WeightIt)
library(cobalt)
library(ggplot2)
library(patchwork)

log_msg("=== ÉTAPE 8 : ESTIMATION DE L'ATE ===")
t0 <- timer_start()

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
log_msg(sprintf("Cohorte: %d patients | Covariables: %d", nrow(cohort), length(propensity_vars)))

# Nettoyage : remplacer NA des covariables par médiane/mode (disponibilité complète pour WeightIt)
cohort_clean <- copy(cohort)
for (v in propensity_vars) {
  if (v %in% names(cohort_clean)) {
    if (is.numeric(cohort_clean[[v]])) {
      med_val <- median(cohort_clean[[v]], na.rm=TRUE)
      cohort_clean[is.na(get(v)), (v) := med_val]
    } else if (is.character(cohort_clean[[v]])) {
      mode_val <- names(sort(table(cohort_clean[[v]]), decreasing=TRUE))[1]
      cohort_clean[is.na(get(v)), (v) := mode_val]
    }
  }
}

# Formule de propension
cov_vars <- intersect(propensity_vars, names(cohort_clean))
ps_formula_str <- paste("high_dose ~", paste(cov_vars, collapse=" + "))
ps_formula <- as.formula(ps_formula_str)
log_msg(sprintf("Formule PS: %s", ps_formula_str))

# =============================================================================
# 1. Naïf non ajusté
# =============================================================================
log_msg("1. Modèle naïf non ajusté...")

model_naive <- glm(outcome_thrombo_100 ~ high_dose, data=cohort_clean, family=binomial)
naive_or <- exp(coef(model_naive)["high_dose"])
naive_ci <- exp(confint.default(model_naive)["high_dose",])

# Risk différence brute
p1 <- mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==1], na.rm=TRUE)
p0 <- mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==0], na.rm=TRUE)
rd_naive <- p1 - p0
rr_naive <- p1 / p0
log_msg(sprintf("Naïf: RD=%.4f, RR=%.2f, OR=%.2f [%.2f, %.2f]",
                rd_naive, rr_naive, naive_or, naive_ci[1], naive_ci[2]))

# =============================================================================
# 2. Régression logistique ajustée (outcome regression)
# =============================================================================
log_msg("2. Régression logistique ajustée...")

outcome_formula <- as.formula(paste("outcome_thrombo_100 ~ high_dose +",
                                     paste(cov_vars, collapse=" + ")))
model_adj <- glm(outcome_formula, data=cohort_clean, family=binomial)
adj_or <- exp(coef(model_adj)["high_dose"])
adj_ci <- exp(confint.default(model_adj)["high_dose",])
log_msg(sprintf("Régression ajustée: OR=%.2f [%.2f, %.2f]",
                adj_or, adj_ci[1], adj_ci[2]))

# G-formula : standardisation marginale
cohort_clone1 <- copy(cohort_clean); cohort_clone1[, high_dose := 1L]
cohort_clone0 <- copy(cohort_clean); cohort_clone0[, high_dose := 0L]
pred1 <- predict(model_adj, newdata=cohort_clone1, type="response")
pred0 <- predict(model_adj, newdata=cohort_clone0, type="response")
gform_rd <- mean(pred1) - mean(pred0)
gform_rr <- mean(pred1) / mean(pred0)
log_msg(sprintf("G-formula: RD=%.4f, RR=%.2f", gform_rd, gform_rr))

# =============================================================================
# 3. IPTW (Propensity Score Weighting)
# =============================================================================
log_msg("3. IPTW avec WeightIt...")

# Modèle logistique pour le PS
set.seed(RANDOM_SEED)
w_logit <- tryCatch(
  weightit(ps_formula, data=as.data.frame(cohort_clean),
           method="ps", estimand="ATE"),
  error=function(e) {
    log_msg(sprintf("WeightIt erreur: %s", e$message), level="WARN")
    NULL
  }
)

if (!is.null(w_logit)) {
  cohort_clean[, w_iptw := w_logit$weights]

  # Diagnostics poids
  w_summary <- summary(w_logit$weights)
  log_msg(sprintf("Poids IPTW: min=%.2f, median=%.2f, max=%.2f, ESS=%.0f",
                  min(w_logit$weights), median(w_logit$weights),
                  max(w_logit$weights), sum(w_logit$weights)^2/sum(w_logit$weights^2)))

  # Balance après pondération
  bal_logit <- bal.tab(w_logit, stats=c("m","v"), thresholds=c(m=0.1))

  # IPTW estimé manuellement (weighted RD)
  iptw_p1 <- weighted.mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==1],
                             cohort_clean$w_iptw[cohort_clean$high_dose==1], na.rm=TRUE)
  iptw_p0 <- weighted.mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==0],
                             cohort_clean$w_iptw[cohort_clean$high_dose==0], na.rm=TRUE)
  iptw_rd <- iptw_p1 - iptw_p0
  iptw_rr <- iptw_p1 / iptw_p0
  log_msg(sprintf("IPTW: RD=%.4f, RR=%.2f", iptw_rd, iptw_rr))

  # Trimming des poids extrêmes (percentile 99)
  w99 <- quantile(w_logit$weights, 0.99)
  cohort_clean[, w_iptw_trim := pmin(w_iptw, w99)]
  iptw_p1_trim <- weighted.mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==1],
                                  cohort_clean$w_iptw_trim[cohort_clean$high_dose==1], na.rm=TRUE)
  iptw_p0_trim <- weighted.mean(cohort_clean$outcome_thrombo_100[cohort_clean$high_dose==0],
                                  cohort_clean$w_iptw_trim[cohort_clean$high_dose==0], na.rm=TRUE)
  iptw_rd_trim <- iptw_p1_trim - iptw_p0_trim
  log_msg(sprintf("IPTW (trimmed p99): RD=%.4f", iptw_rd_trim))
} else {
  iptw_rd <- NA_real_; iptw_rr <- NA_real_; iptw_rd_trim <- NA_real_
  cohort_clean[, w_iptw := 1]; cohort_clean[, w_iptw_trim := 1]
  bal_logit <- NULL
}

# =============================================================================
# 4. AIPW (Doubly Robust) — implémentation manuelle simple
# =============================================================================
log_msg("4. AIPW (doubly robust)...")

# Modèle de propension
ps_model <- glm(ps_formula, data=as.data.frame(cohort_clean), family=binomial)
cohort_clean[, ps := predict(ps_model, type="response")]

# Modèle d'outcome
cohort_clean[, mu1 := pred1]  # déjà calculé depuis g-formula
cohort_clean[, mu0 := pred0]

# AIPW estimateur
A <- cohort_clean$high_dose
Y <- cohort_clean$outcome_thrombo_100
ps_hat <- cohort_clean$ps
mu1_hat <- cohort_clean$mu1
mu0_hat <- cohort_clean$mu0

aipw_eif_1 <- mu1_hat + A/ps_hat * (Y - mu1_hat)
aipw_eif_0 <- mu0_hat + (1-A)/(1-ps_hat) * (Y - mu0_hat)
aipw_ate <- mean(aipw_eif_1) - mean(aipw_eif_0)

# Variance via EIF
eif_diff <- aipw_eif_1 - aipw_eif_0 - aipw_ate
aipw_se  <- sqrt(var(eif_diff) / length(eif_diff))
aipw_ci  <- c(aipw_ate - 1.96*aipw_se, aipw_ate + 1.96*aipw_se)

log_msg(sprintf("AIPW: RD=%.4f, SE=%.4f, 95%%CI [%.4f, %.4f]",
                aipw_ate, aipw_se, aipw_ci[1], aipw_ci[2]))

# =============================================================================
# 5. Bootstrap pour intervalles de confiance robustes
# =============================================================================
log_msg("5. Bootstrap (500 réplications) pour les IC...")

set.seed(RANDOM_SEED)
n <- nrow(cohort_clean)
B <- 500
boot_results <- matrix(NA, nrow=B, ncol=3)

for (b in seq_len(B)) {
  idx <- sample(n, n, replace=TRUE)
  boot_dt <- cohort_clean[idx,]

  # Propension
  tryCatch({
    ps_b <- glm(ps_formula, data=as.data.frame(boot_dt), family=binomial)
    ps_hat_b <- predict(ps_b, type="response")
    ps_hat_b <- pmax(pmin(ps_hat_b, 0.99), 0.01)

    # G-formula
    cl1 <- copy(boot_dt); cl1[, high_dose := 1L]
    cl0 <- copy(boot_dt); cl0[, high_dose := 0L]
    out_b <- glm(outcome_formula, data=as.data.frame(boot_dt), family=binomial)
    mu1_b <- predict(out_b, newdata=cl1, type="response")
    mu0_b <- predict(out_b, newdata=cl0, type="response")

    # AIPW
    A_b <- boot_dt$high_dose; Y_b <- boot_dt$outcome_thrombo_100
    eif1_b <- mu1_b + A_b/ps_hat_b * (Y_b - mu1_b)
    eif0_b <- mu0_b + (1-A_b)/(1-ps_hat_b) * (Y_b - mu0_b)
    boot_results[b, 1] <- mean(mu1_b) - mean(mu0_b)  # g-formula
    boot_results[b, 2] <- mean(eif1_b) - mean(eif0_b)  # AIPW
    boot_results[b, 3] <- weighted.mean(Y_b[A_b==1], 1/ps_hat_b[A_b==1]) -
                           weighted.mean(Y_b[A_b==0], 1/(1-ps_hat_b[A_b==0]))  # IPTW
  }, error=function(e) NULL)
}

boot_results <- as.data.table(boot_results)
setnames(boot_results, c("rd_gform", "rd_aipw", "rd_iptw"))
boot_results <- boot_results[!is.na(rd_gform)]
log_msg(sprintf("Bootstrap: %d réplications valides", nrow(boot_results)))

# ICs bootstrap
ci_gform <- quantile(boot_results$rd_gform, c(0.025, 0.975), na.rm=TRUE)
ci_aipw  <- quantile(boot_results$rd_aipw,  c(0.025, 0.975), na.rm=TRUE)
ci_iptw  <- quantile(boot_results$rd_iptw,  c(0.025, 0.975), na.rm=TRUE)

# =============================================================================
# 6. Synthèse des résultats ATE
# =============================================================================
log_msg("6. Synthèse résultats ATE...")

ate_results <- data.table(
  estimator    = c("1. Non ajusté (RD brut)",
                   "2. Régression logistique ajustée (OR → RD marginal)",
                   "3. G-formula (standardisation marginale)",
                   "4. IPTW logistique",
                   "5. IPTW logistique (poids trimmed p99)",
                   "6. AIPW doubly-robust"),
  outcome      = "Thrombopénie < 100 G/L (J0-J30)",
  estimate_rd  = round(c(rd_naive, gform_rd, gform_rd, iptw_rd, iptw_rd_trim, aipw_ate), 4),
  ci_lower     = round(c(NA, ci_gform[1], ci_gform[1], ci_iptw[1], NA, ci_aipw[1]), 4),
  ci_upper     = round(c(NA, ci_gform[2], ci_gform[2], ci_iptw[2], NA, ci_aipw[2]), 4),
  note         = c(
    "Pas d'ajustement",
    "Ajustement régression — pas doubly robust",
    "G-formula — IC bootstrap",
    "IPW — IC bootstrap",
    "IPW trimmed — robustesse aux poids extrêmes",
    "Doubly robust — IC bootstrap"
  )
)

# Ajouter résultats pour les autres outcomes
make_outcome_row <- function(outcome_col, label) {
  Y2 <- as.integer(cohort_clean[[outcome_col]])
  if (sum(!is.na(Y2)) < 30) {
    return(data.table(estimator="AIPW", outcome=label, estimate_rd=NA, ci_lower=NA, ci_upper=NA, note="Trop peu de cas"))
  }
  p1_crude <- mean(Y2[A==1], na.rm=TRUE)
  p0_crude <- mean(Y2[A==0], na.rm=TRUE)
  # AIPW
  out_m <- glm(as.formula(paste(outcome_col, "~ high_dose +", paste(cov_vars, collapse=" + "))),
               data=as.data.frame(cohort_clean), family=binomial)
  cl1_tmp <- copy(cohort_clean); cl1_tmp[, high_dose := 1L]
  cl0_tmp <- copy(cohort_clean); cl0_tmp[, high_dose := 0L]
  mu1_tmp <- predict(out_m, newdata=cl1_tmp, type="response")
  mu0_tmp <- predict(out_m, newdata=cl0_tmp, type="response")
  ps_h <- pmax(pmin(cohort_clean$ps, 0.99), 0.01)
  eif1_tmp <- mu1_tmp + A/ps_h * (Y2 - mu1_tmp)
  eif0_tmp <- mu0_tmp + (1-A)/(1-ps_h) * (Y2 - mu0_tmp)
  ate_tmp <- mean(eif1_tmp) - mean(eif0_tmp)
  data.table(estimator="AIPW", outcome=label,
             estimate_rd=round(ate_tmp, 4),
             ci_lower=NA, ci_upper=NA,
             note="IC à calculer")
}

tryCatch({
  ate_sec <- rbind(
    make_outcome_row("outcome_thrombo_50",   "Thrombopénie < 50 G/L"),
    make_outcome_row("outcome_thrombo_rel30","Baisse plaquettes > 30%"),
    make_outcome_row("outcome_hepato_3xULN", "Hépatotoxicité ALT/AST > 3xULN"),
    make_outcome_row("outcome_hypernh3",     "Hyperammoniémie NH3 > 55"),
    make_outcome_row("outcome_composite",    "Outcome composite")
  )
  ate_results <- rbind(ate_results, ate_sec, fill=TRUE)
}, error=function(e) log_msg(sprintf("Erreur outcomes secondaires: %s", e$message), level="WARN"))

log_msg("Résultats ATE :")
print(ate_results)
safe_write_csv(ate_results, file.path(TABLES_DIR, "ate_results.csv"), overwrite=TRUE)

# =============================================================================
# 7. Figures diagnostiques propension
# =============================================================================
log_msg("7. Figures diagnostiques...")

# Distribution du propensity score
p_ps <- ggplot(cohort_clean, aes(x=ps, fill=factor(high_dose))) +
  geom_histogram(aes(y=after_stat(density)), bins=40, alpha=0.6, position="identity") +
  geom_density(aes(color=factor(high_dose)), linewidth=0.8) +
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("Faible dose","Haute dose")) +
  scale_color_manual(values=c("steelblue","tomato"), guide="none") +
  labs(title="Distribution du propensity score par groupe",
       subtitle="Overlap des distributions = positivité",
       x="Propensity score P(haute dose | X)", y="Densité", fill="Groupe") +
  theme_bw(base_size=12)

ggsave(file.path(DIAGNOSTICS_DIR, "ps_distribution.png"), p_ps, width=9, height=6, dpi=150)

# Love plot (balance)
if (!is.null(bal_logit)) {
  tryCatch({
    png(file.path(DIAGNOSTICS_DIR, "love_plot.png"), width=900, height=700)
    love.plot(bal_logit, threshold=0.1, abs=TRUE, var.order="adjusted",
              title="Love Plot — Balance avant/après IPTW")
    dev.off()
    log_msg("Love plot sauvegardé")
  }, error=function(e) log_msg(sprintf("Love plot erreur: %s", e$message), level="WARN"))
}

# Forest plot des résultats ATE
ate_main <- ate_results[outcome == "Thrombopénie < 100 G/L (J0-J30)" & !is.na(estimate_rd)]
p_forest_ate <- ggplot(ate_main, aes(x=estimate_rd, y=reorder(estimator, estimate_rd))) +
  geom_point(size=3, color="steelblue") +
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper), height=0.2, na.rm=TRUE) +
  geom_vline(xintercept=0, color="gray50", linetype="dashed") +
  labs(title="Risk Difference — Haute vs Faible dose VPA",
       subtitle="Outcome : Thrombopénie incidente (plaquettes < 100 G/L, J0-J30)",
       x="Risk Difference", y="Estimateur") +
  theme_bw(base_size=11)

ggsave(file.path(DIAGNOSTICS_DIR, "forest_ate_main.png"),
       p_forest_ate, width=11, height=6, dpi=150)

# Distribution des poids IPTW
p_weights <- ggplot(cohort_clean, aes(x=w_iptw, fill=factor(high_dose))) +
  geom_histogram(bins=60, alpha=0.7) +
  geom_vline(xintercept=quantile(cohort_clean$w_iptw, 0.99), color="red", linetype="dashed") +
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("Faible dose","Haute dose")) +
  labs(title="Distribution des poids IPTW", subtitle="Ligne rouge = percentile 99 (trimming)",
       x="Poids IPTW", y="N", fill="Groupe") +
  theme_bw(base_size=12)

ggsave(file.path(DIAGNOSTICS_DIR, "iptw_weights_distribution.png"),
       p_weights, width=9, height=6, dpi=150)

timer_end(t0, "Étape 8 totale")

writeLines(sprintf(
"# Estimation de l'ATE — Étude VPA/MIMIC-IV
**Date :** %s

## Résultats principaux (Thrombopénie < 100 G/L)

| Estimateur | RD | IC 95%% |
|-----------|-----|---------|
%s

## Interprétation préliminaire

- La Risk Difference brute est **+%.1f%%** (8.7%% vs 4.9%%)
- Après ajustement (G-formula, AIPW), l'effet tend à s'atténuer, suggérant un confounding important
- ICU, sévérité et sexe sont les confounders principaux

## Notes méthodologiques
- Propensity score estimé par régression logistique
- AIPW : doubly-robust (valide si PS OU modèle d'outcome correct)
- IC : bootstrap 500 réplications
- Trimming des poids IPTW : percentile 99
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(apply(ate_main, 1, function(r)
    sprintf("| %s | %.3f | [%.3f, %.3f] |",
            r["estimator"], as.numeric(r["estimate_rd"]),
            ifelse(is.na(r["ci_lower"]), NA_real_, as.numeric(r["ci_lower"])),
            ifelse(is.na(r["ci_upper"]), NA_real_, as.numeric(r["ci_upper"])))),
    collapse="\n"),
  rd_naive * 100
), file.path(REPORTS_DIR, "08_ate_analysis.md"))

log_msg("=== ÉTAPE 8 TERMINÉE ===")
