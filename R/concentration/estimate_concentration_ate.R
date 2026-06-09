# =============================================================================
# estimate_concentration_ate.R
# Conc-3 — ATE estimation using plasma VPA concentrations
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

log_msg("=== ÉTAPE 8B.6 : ATE CONCENTRATION-BASED ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR, DIAGNOSTICS_DIR)

cohort_conc <- read_parquet_dt(file.path(DATA_FINAL, "analytic_cohort_concentration.parquet"))
log_msg(sprintf("Cohorte concentration: %d patients | Haute conc: %d (%.1f%%)",
                nrow(cohort_conc), sum(cohort_conc$high_conc_75),
                mean(cohort_conc$high_conc_75)*100))

# =============================================================================
# Variables du propensity score pour l'analyse concentration-based
# (= déterminants de la concentration élevée, différent du PS dose-based)
# Inclure : déterminants dose ET déterminants PK/TDM
# =============================================================================
ps_vars_conc <- c(
  # Dose prescrite (principal déterminant de la concentration)
  "daily_dose_main",
  # Facteurs pharmacocinétiques
  "age", "female", "weight_kg_baseline", "alb_baseline", "crea_baseline",
  # Contexte clinique
  "in_icu", "sepsis", "emergency",
  # Indications
  "epilepsy_dx", "bipolar_dx",
  # Route d'administration
  "vpa_iv",
  # Comédications
  "comed_topiramate", "comed_carbapenem", "comed_heparin",
  "comed_steroid", "comed_antipsy", "comed_other_aed",
  # Baseline labs
  "plt_baseline", "alt_baseline", "ast_baseline", "liver_baseline_elevated",
  # Autres
  "pre_hosp_days_winsorized"
)
ps_vars_conc <- intersect(ps_vars_conc, names(cohort_conc))

# Imputation simple des NA par médiane/mode
cohort_clean <- copy(cohort_conc)
for (v in ps_vars_conc) {
  if (is.numeric(cohort_clean[[v]])) {
    med_v <- median(cohort_clean[[v]], na.rm=TRUE)
    cohort_clean[is.na(get(v)), (v) := med_v]
  } else {
    mode_v <- names(sort(table(cohort_clean[[v]]), decreasing=TRUE))[1]
    cohort_clean[is.na(get(v)), (v) := mode_v]
  }
}

# Formule propension
ps_formula_str <- paste("high_conc_75 ~", paste(ps_vars_conc, collapse=" + "))
ps_formula <- as.formula(ps_formula_str)
log_msg(sprintf("PS formula (%d variables): %s", length(ps_vars_conc), ps_formula_str))

A <- cohort_clean$high_conc_75
Y <- as.integer(cohort_clean$outcome_thrombo_100)
n <- nrow(cohort_clean)

# =============================================================================
# 1. Naïf non ajusté
# =============================================================================
log_msg("1. Modèle naïf...")
p1_naive <- mean(Y[A==1], na.rm=TRUE)
p0_naive <- mean(Y[A==0], na.rm=TRUE)
rd_naive <- p1_naive - p0_naive
rr_naive <- p1_naive / p0_naive
log_msg(sprintf("Naïf: P(Y=1|haute) = %.3f, P(Y=1|basse) = %.3f, RD = %.4f, RR = %.2f",
                p1_naive, p0_naive, rd_naive, rr_naive))

# =============================================================================
# 2. G-formula (outcome regression)
# =============================================================================
log_msg("2. G-formula (outcome regression)...")
outcome_formula <- as.formula(paste("outcome_thrombo_100 ~ high_conc_75 +",
                                     paste(ps_vars_conc, collapse=" + ")))
model_out <- glm(outcome_formula, data=as.data.frame(cohort_clean), family=binomial)

cl1 <- copy(cohort_clean); cl1[, high_conc_75 := 1L]
cl0 <- copy(cohort_clean); cl0[, high_conc_75 := 0L]
mu1 <- predict(model_out, newdata=cl1, type="response")
mu0 <- predict(model_out, newdata=cl0, type="response")
gform_rd <- mean(mu1) - mean(mu0)
gform_rr <- mean(mu1) / mean(mu0)
log_msg(sprintf("G-formula: RD=%.4f, RR=%.2f", gform_rd, gform_rr))

# =============================================================================
# 3. IPTW
# =============================================================================
log_msg("3. IPTW avec WeightIt...")
set.seed(RANDOM_SEED)
w_logit <- tryCatch(
  weightit(ps_formula, data=as.data.frame(cohort_clean), method="ps", estimand="ATE"),
  error=function(e) { log_msg(sprintf("WARN: %s", e$message)); NULL }
)

if (!is.null(w_logit)) {
  cohort_clean[, w_iptw := w_logit$weights]
  w_ess <- sum(w_logit$weights)^2 / sum(w_logit$weights^2)
  log_msg(sprintf("IPTW poids: min=%.2f, median=%.2f, max=%.2f, ESS=%.0f/%.0f",
                  min(w_logit$weights), median(w_logit$weights),
                  max(w_logit$weights), w_ess, n))

  iptw_p1 <- weighted.mean(Y[A==1], w_logit$weights[A==1])
  iptw_p0 <- weighted.mean(Y[A==0], w_logit$weights[A==0])
  iptw_rd <- iptw_p1 - iptw_p0
  iptw_rr <- iptw_p1 / iptw_p0
  log_msg(sprintf("IPTW: RD=%.4f, RR=%.2f", iptw_rd, iptw_rr))

  w99 <- quantile(w_logit$weights, 0.99)
  cohort_clean[, w_iptw_trim := pmin(w_iptw, w99)]
  iptw_p1_t <- weighted.mean(Y[A==1], cohort_clean$w_iptw_trim[A==1])
  iptw_p0_t <- weighted.mean(Y[A==0], cohort_clean$w_iptw_trim[A==0])
  iptw_rd_trim <- iptw_p1_t - iptw_p0_t
  log_msg(sprintf("IPTW (trim p99): RD=%.4f", iptw_rd_trim))

  # Balance diagnostics
  bal_conc <- bal.tab(w_logit, stats=c("m","v"), thresholds=c(m=0.1))
} else {
  cohort_clean[, w_iptw := 1]; cohort_clean[, w_iptw_trim := 1]
  iptw_rd <- NA; iptw_rr <- NA; iptw_rd_trim <- NA
  w_ess <- NA; w_logit <- NULL; bal_conc <- NULL
}

# =============================================================================
# 4. AIPW (doubly robust)
# =============================================================================
log_msg("4. AIPW (doubly robust)...")
ps_model <- glm(ps_formula, data=as.data.frame(cohort_clean), family=binomial)
ps_hat <- predict(ps_model, type="response")
ps_hat <- pmax(pmin(ps_hat, 0.99), 0.01)

aipw_eif_1 <- mu1 + A / ps_hat * (Y - mu1)
aipw_eif_0 <- mu0 + (1-A) / (1-ps_hat) * (Y - mu0)
aipw_ate <- mean(aipw_eif_1) - mean(aipw_eif_0)

eif_diff  <- aipw_eif_1 - aipw_eif_0 - aipw_ate
aipw_se   <- sqrt(var(eif_diff) / length(eif_diff))
aipw_ci   <- c(aipw_ate - 1.96*aipw_se, aipw_ate + 1.96*aipw_se)

log_msg(sprintf("AIPW: RD=%.4f, SE=%.4f, 95%%CI [%.4f, %.4f]",
                aipw_ate, aipw_se, aipw_ci[1], aipw_ci[2]))

# =============================================================================
# 5. Bootstrap pour IC
# =============================================================================
log_msg("5. Bootstrap (500 réplications)...")
set.seed(RANDOM_SEED)
B <- 500
boot_mat <- matrix(NA, B, 3)

for (b in seq_len(B)) {
  idx <- sample(n, n, replace=TRUE)
  bd  <- cohort_clean[idx,]
  tryCatch({
    ps_b <- glm(ps_formula, data=as.data.frame(bd), family=binomial)
    ps_b_hat <- pmax(pmin(predict(ps_b, type="response"), 0.99), 0.01)
    out_b <- glm(outcome_formula, data=as.data.frame(bd), family=binomial)
    b1  <- copy(bd); b1[, high_conc_75 := 1L]
    b0  <- copy(bd); b0[, high_conc_75 := 0L]
    m1b <- predict(out_b, newdata=b1, type="response")
    m0b <- predict(out_b, newdata=b0, type="response")
    Ab  <- bd$high_conc_75; Yb <- as.integer(bd$outcome_thrombo_100)
    e1b <- m1b + Ab/ps_b_hat * (Yb - m1b)
    e0b <- m0b + (1-Ab)/(1-ps_b_hat) * (Yb - m0b)
    boot_mat[b, 1] <- mean(m1b) - mean(m0b)
    boot_mat[b, 2] <- mean(e1b) - mean(e0b)
    boot_mat[b, 3] <- weighted.mean(Yb[Ab==1], 1/ps_b_hat[Ab==1]) -
                      weighted.mean(Yb[Ab==0], 1/(1-ps_b_hat[Ab==0]))
  }, error=function(e) NULL)
}
boot_dt <- as.data.table(boot_mat)
setnames(boot_dt, c("rd_gform","rd_aipw","rd_iptw"))
boot_dt <- boot_dt[!is.na(rd_gform)]
log_msg(sprintf("Bootstrap valide: %d/%d réplications", nrow(boot_dt), B))

ci_gf  <- quantile(boot_dt$rd_gform, c(0.025,0.975), na.rm=TRUE)
ci_ai  <- quantile(boot_dt$rd_aipw,  c(0.025,0.975), na.rm=TRUE)
ci_ip  <- quantile(boot_dt$rd_iptw,  c(0.025,0.975), na.rm=TRUE)

# =============================================================================
# 6. Synthèse des estimateurs — outcome principal (thrombopénie)
# =============================================================================
ate_main <- data.table(
  estimator   = c("1. Non ajusté (RD brut)",
                  "2. G-formula (standardisation marginale)",
                  "3. IPTW logistique",
                  "4. IPTW (poids trimmed p99)",
                  "5. AIPW doubly-robust"),
  outcome     = "Thrombopénie < 100 G/L",
  estimate_rd = round(c(rd_naive, gform_rd, iptw_rd, iptw_rd_trim, aipw_ate), 4),
  ci_lower    = round(c(NA, ci_gf[1], ci_ip[1], NA, ci_ai[1]), 4),
  ci_upper    = round(c(NA, ci_gf[2], ci_ip[2], NA, ci_ai[2]), 4),
  p_exposed   = round(c(p1_naive, mean(mu1), iptw_p1, NA, NA), 4),
  p_unexposed = round(c(p0_naive, mean(mu0), iptw_p0, NA, NA), 4)
)

# Outcomes secondaires via AIPW rapide
make_aipw_row <- function(oc, label) {
  Y2 <- as.integer(cohort_clean[[oc]])
  if (sum(!is.na(Y2)) < 20 || sum(Y2, na.rm=TRUE) < 10) {
    return(data.table(estimator="AIPW", outcome=label, estimate_rd=NA, ci_lower=NA, ci_upper=NA,
                      p_exposed=NA, p_unexposed=NA, note="Trop peu de cas"))
  }
  tryCatch({
    out_m <- glm(as.formula(paste(oc, "~ high_conc_75 +", paste(ps_vars_conc, collapse=" + "))),
                 data=as.data.frame(cohort_clean), family=binomial)
    cl1_t <- copy(cohort_clean); cl1_t[, high_conc_75 := 1L]
    cl0_t <- copy(cohort_clean); cl0_t[, high_conc_75 := 0L]
    m1_t <- predict(out_m, newdata=cl1_t, type="response")
    m0_t <- predict(out_m, newdata=cl0_t, type="response")
    e1_t <- m1_t + A/ps_hat * (Y2 - m1_t)
    e0_t <- m0_t + (1-A)/(1-ps_hat) * (Y2 - m0_t)
    ate_t <- mean(e1_t) - mean(e0_t)
    # Bootstrap rapide (200 reps)
    set.seed(RANDOM_SEED + 1)
    boot2 <- sapply(1:200, function(b) {
      tryCatch({
        idx2 <- sample(n, n, replace=TRUE)
        bd2  <- cohort_clean[idx2,]
        ps2  <- pmax(pmin(predict(glm(ps_formula, data=as.data.frame(bd2), family=binomial),
                                   type="response"), 0.99), 0.01)
        om2  <- glm(as.formula(paste(oc, "~ high_conc_75 +", paste(ps_vars_conc, collapse=" + "))),
                    data=as.data.frame(bd2), family=binomial)
        c1 <- copy(bd2); c1[, high_conc_75:=1L]
        c0 <- copy(bd2); c0[, high_conc_75:=0L]
        m1b2 <- predict(om2, newdata=c1, type="response")
        m0b2 <- predict(om2, newdata=c0, type="response")
        Ab2 <- bd2$high_conc_75; Yb2 <- as.integer(bd2[[oc]])
        mean(m1b2 + Ab2/ps2*(Yb2-m1b2)) - mean(m0b2 + (1-Ab2)/(1-ps2)*(Yb2-m0b2))
      }, error=function(e) NA_real_)
    })
    ci_t <- quantile(boot2, c(0.025,0.975), na.rm=TRUE)
    data.table(estimator="AIPW", outcome=label, estimate_rd=round(ate_t,4),
               ci_lower=round(ci_t[1],4), ci_upper=round(ci_t[2],4),
               p_exposed=round(mean(m1_t),4), p_unexposed=round(mean(m0_t),4), note="IC bootstrap 200 reps")
  }, error=function(e) {
    data.table(estimator="AIPW", outcome=label, estimate_rd=NA, ci_lower=NA, ci_upper=NA,
               p_exposed=NA, p_unexposed=NA, note=as.character(e$message))
  })
}

ate_sec <- rbindlist(list(
  make_aipw_row("outcome_thrombo_50",   "Thrombopénie < 50 G/L"),
  make_aipw_row("outcome_hepato_3xULN", "Hépatotoxicité ALT/AST > 3xULN"),
  make_aipw_row("outcome_hypernh3",     "Hyperammoniémie NH3 > 55"),
  make_aipw_row("outcome_composite",    "Outcome composite")
), fill=TRUE)

ate_all <- rbind(ate_main, ate_sec, fill=TRUE)
log_msg("Résultats ATE concentration-based:")
print(ate_all[, .(estimator, outcome, estimate_rd, ci_lower, ci_upper)])
safe_write_csv(ate_all, file.path(TABLES_DIR, "ate_results_concentration.csv"), overwrite=TRUE)

# =============================================================================
# 7. Diagnostics propension
# =============================================================================
log_msg("7. Diagnostics propension...")

cohort_clean[, ps_conc := ps_hat]

# Distribution du PS par groupe
p_ps <- ggplot(cohort_clean, aes(x=ps_conc, fill=factor(high_conc_75))) +
  geom_histogram(aes(y=after_stat(density)), bins=40, alpha=0.6, position="identity") +
  geom_density(aes(color=factor(high_conc_75)), linewidth=0.8) +
  scale_fill_manual(values=c("steelblue","tomato"),
                    labels=c("Basse conc. (<75)","Haute conc. (≥75)")) +
  scale_color_manual(values=c("steelblue","tomato"), guide="none") +
  labs(title="Distribution du propensity score (concentration-based)",
       subtitle=sprintf("N=%d | ESS=%.0f/%.0f", n, w_ess, n),
       x="P(haute concentration ≥ 75 µg/mL | X)", y="Densité", fill="Groupe") +
  theme_bw(base_size=12)
ggsave(file.path(DIAGNOSTICS_DIR, "propensity_diagnostics_concentration.png"),
       p_ps, width=9, height=6, dpi=150)

# Forest plot des estimateurs
ate_plot <- ate_all[!is.na(estimate_rd) & estimator %in%
  c("1. Non ajusté (RD brut)","2. G-formula (standardisation marginale)",
    "3. IPTW logistique","5. AIPW doubly-robust")]
ate_plot_prim <- ate_plot[outcome == "Thrombopénie < 100 G/L"]

p_ate <- ggplot(ate_plot_prim, aes(x=estimate_rd*100, y=reorder(estimator, -estimate_rd))) +
  geom_point(size=3, color="tomato") +
  geom_errorbarh(aes(xmin=ci_lower*100, xmax=ci_upper*100), height=0.3,
                  color="tomato", linewidth=0.8) +
  geom_vline(xintercept=0, linetype="dashed", color="gray30") +
  scale_x_continuous(labels=function(x) paste0(x, "%")) +
  labs(title="Estimateurs ATE — Thrombopénie (concentration-based)",
       subtitle="Exposition: première concentration VPA ≥ 75 µg/mL vs < 75 µg/mL",
       x="Risk Difference (percentage points)", y="Estimateur") +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "ate_forest_concentration.png"),
       p_ate, width=10, height=5, dpi=150)

# =============================================================================
# Sauvegarde PS pour usage aval
# =============================================================================
arrow::write_parquet(as.data.frame(cohort_clean[, .(subject_id, ps_conc, w_iptw, w_iptw_trim)]),
                      file.path(DATA_INTERMEDIATE, "ps_concentration.parquet"))

# =============================================================================
# Rapport 8B.6
# =============================================================================
aipw_row_prim <- ate_all[estimator == "5. AIPW doubly-robust" &
                            outcome == "Thrombopénie < 100 G/L"]

report_8b6 <- sprintf(
'# 8B.6 Estimation ATE — Concentration-Based
**Date :** %s

---

## Population

- Sous-cohorte monitorée : **%d patients**
- Haute concentration (≥ 75 µg/mL) : **%d (%.1f%%)**
- Basse concentration (< 75 µg/mL) : **%d (%.1f%%)**

## Variables du propensity score

%d variables : %s

## Diagnostics du propensity score

- ESS (effective sample size) : **%.0f / %d**
- Poids IPTW : médiane = %.2f, max = %.2f
- Stabilisation : trimming au percentile 99 (seuil = %.2f)

## Résultats ATE — Outcome principal : Thrombopénie < 100 G/L

| Estimateur | RD | IC 95%% |
|------------|-----|---------|
| Non ajusté (brut) | %.1f%% | — |
| G-formula | %.1f%% | [%.1f%%, %.1f%%] |
| IPTW logistique | %.1f%% | [%.1f%%, %.1f%%] |
| AIPW (doubly robust) | **%.1f%%** | **[%.1f%%, %.1f%%]** |

## Outcomes secondaires (AIPW)

%s

## Interprétation

L\'estimateur AIPW doubly-robust (résultats primaires) fournit l\'estimation la plus robuste
de l\'association causale entre concentration VPA ≥ 75 µg/mL et thrombopénie incidente.

**Limites critiques :**
1. L\'ESS (%.0f/%.0f = %.0f%%) reflète les déséquilibres entre groupes — pondération imparfaite
2. La sélection des patients monitorés est fortement informative (confounding by indication du TDM)
3. Le timing trough/peak est inconnu — variabilité non contrôlée
4. Résultats à interpréter comme exploratoires

**Bootstrap :** %d réplications valides sur %d
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(cohort_conc),
  sum(cohort_conc$high_conc_75), mean(cohort_conc$high_conc_75)*100,
  sum(cohort_conc$high_conc_75==0), mean(cohort_conc$high_conc_75==0)*100,
  length(ps_vars_conc),
  paste(ps_vars_conc, collapse=", "),
  w_ess, n,
  if(!is.null(w_logit)) median(w_logit$weights) else NA,
  if(!is.null(w_logit)) max(w_logit$weights) else NA,
  if(!is.null(w_logit)) quantile(w_logit$weights, 0.99) else NA,
  rd_naive*100,
  gform_rd*100, ci_gf[1]*100, ci_gf[2]*100,
  iptw_rd*100, ci_ip[1]*100, ci_ip[2]*100,
  aipw_ate*100, aipw_ci[1]*100, aipw_ci[2]*100,
  {
    sec_rows <- ate_sec[!is.na(estimate_rd)]
    if (nrow(sec_rows) == 0) "Aucun outcome secondaire estimable"
    else paste(sprintf("| %s | %.1f%% | [%.1f%%, %.1f%%] |",
                       sec_rows$outcome, sec_rows$estimate_rd*100,
                       sec_rows$ci_lower*100, sec_rows$ci_upper*100),
               collapse="\n")
  },
  w_ess, n, w_ess/n*100,
  nrow(boot_dt), B
)
writeLines(report_8b6, file.path(REPORTS_DIR, "08B_ate_concentration.md"))
log_msg("Rapport 8B.6 écrit")

timer_end(t0, "8B.6 ATE concentration")
log_msg("=== ÉTAPE 8B.6 TERMINÉE ===")
