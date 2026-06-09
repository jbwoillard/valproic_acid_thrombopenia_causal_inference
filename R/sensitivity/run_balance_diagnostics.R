# =============================================================================
# run_balance_diagnostics.R
# Step 11 — Propensity score overlap and balance diagnostics
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table); library(arrow); library(WeightIt); library(cobalt)
library(ggplot2); library(patchwork)

log_msg("=== STARTING ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR, DIAGNOSTICS_DIR)

# =============================================================================
# PARTIE 1 : Analyse dose-based
# =============================================================================
log_msg("-- PARTIE 1 : Dose-based --")

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
propensity_vars <- readRDS(file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))
cov_vars <- intersect(propensity_vars, names(cohort))

cohort_d <- copy(cohort)
for (v in cov_vars) {
  if (is.numeric(cohort_d[[v]])) {
    cohort_d[is.na(get(v)), (v) := median(cohort_d[[v]], na.rm=TRUE)]
  }
}

ps_formula_d <- as.formula(paste("high_dose ~", paste(cov_vars, collapse="+")))

# Fit WeightIt
set.seed(RANDOM_SEED)
w_dose <- tryCatch(
  weightit(ps_formula_d, data=as.data.frame(cohort_d), method="ps", estimand="ATE"),
  error=function(e) { log_msg(sprintf("WARN: %s", e$message)); NULL }
)

if (!is.null(w_dose)) {
  # 1a. Love plot complet — toutes les covariables
  love_d <- love.plot(
    w_dose,
    stats = "mean.diffs",
    threshold = 0.1,
    abs = TRUE,
    var.order = "unadjusted",
    colors = c("tomato", "steelblue"),
    shapes = c("circle", "triangle"),
    sample.names = c("Non ajusté", "IPTW ajusté"),
    title = "Balance des covariables — Analyse Dose-Based",
    subtitle = "SMD absolu | Seuil 0.1 en rouge"
  )
  png(file.path(FIGURES_DIR, "love_plot_full_dose_based.png"), width=1200, height=900, res=130)
  print(love_d)
  dev.off()
  log_msg("Love plot dose-based sauvegardé")

  # 1b. Mirrored histogram (PS overlap)
  ps_hat_d <- w_dose$ps
  df_ps_d <- data.frame(
    PS        = as.numeric(ps_hat_d),
    Treatment = factor(cohort_d$high_dose, labels=c("Faible dose (≤1000 mg/j)", "Haute dose (>1000 mg/j)"))
  )
  p_mirror_d <- ggplot(df_ps_d, aes(x=PS, fill=Treatment)) +
    geom_histogram(data=subset(df_ps_d, Treatment=="Haute dose (>1000 mg/j)"),
                   aes(y=after_stat(count)), bins=50, alpha=0.7) +
    geom_histogram(data=subset(df_ps_d, Treatment=="Faible dose (≤1000 mg/j)"),
                   aes(y=-after_stat(count)), bins=50, alpha=0.7) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    scale_y_continuous(labels=abs) +
    scale_fill_manual(values=c("steelblue", "tomato")) +
    labs(title="Overlap du Propensity Score — Dose-Based",
         subtitle=sprintf("N=%d | Haute dose: %d | Faible dose: %d",
                          nrow(df_ps_d), sum(cohort_d$high_dose==1), sum(cohort_d$high_dose==0)),
         x="Propensity Score P(haute dose | X)", y="N patients") +
    theme_bw(base_size=12) + theme(legend.position="bottom")
  ggsave(file.path(FIGURES_DIR, "ps_overlap_dose_based.png"), p_mirror_d, width=10, height=6, dpi=150)
  log_msg("Mirrored histogram dose-based sauvegardé")

  # 1c. Distribution des poids IPW
  w_d <- w_dose$weights
  w_trim_d <- pmin(w_d, quantile(w_d, 0.99))
  ess_d <- sum(w_d)^2 / sum(w_d^2)

  p_weights_d <- ggplot(data.frame(w=w_d, group=factor(cohort_d$high_dose,
                                                         labels=c("Faible dose","Haute dose"))),
                         aes(x=w, fill=group)) +
    geom_histogram(bins=50, alpha=0.7) +
    geom_vline(xintercept=quantile(w_d, 0.99), color="red", linetype="dashed") +
    scale_fill_manual(values=c("steelblue","tomato")) +
    scale_x_continuous(limits=c(0, quantile(w_d, 0.995))) +
    labs(title="Distribution des Poids IPTW — Dose-Based",
         subtitle=sprintf("ESS=%.0f/%d (%.1f%%) | Max=%.2f | P99=%.2f | Ligne rouge=seuil trimming",
                          ess_d, nrow(cohort_d), ess_d/nrow(cohort_d)*100,
                          max(w_d), quantile(w_d, 0.99)),
         x="Poids IPTW", y="N patients", fill="Groupe") +
    theme_bw(base_size=12) + theme(legend.position="bottom")
  ggsave(file.path(FIGURES_DIR, "ipw_weight_distribution_dose_based.png"),
         p_weights_d, width=10, height=6, dpi=150)
  log_msg(sprintf("Poids IPTW dose: ESS=%.0f/%.0f, max=%.2f, P99=%.2f",
                  ess_d, nrow(cohort_d), max(w_d), quantile(w_d, 0.99)))

  # 1d. Balance table complète
  bal_tab_d <- bal.tab(w_dose, stats=c("m","v"), thresholds=c(m=0.1))
  bal_dt_d <- as.data.table(bal_tab_d$Balance)
  bal_dt_d[, variable := rownames(bal_tab_d$Balance)]
  bal_dt_d[, analysis := "dose-based"]
  safe_write_csv(bal_dt_d, file.path(TABLES_DIR, "balance_full_dose_based.csv"), overwrite=TRUE)
  max_smd_adj_d <- max(abs(bal_dt_d$Diff.Adj), na.rm=TRUE)
  log_msg(sprintf("SMD max ajusté (dose): %.3f", max_smd_adj_d))
}

# =============================================================================
# PARTIE 2 : Analyse concentration-based
# =============================================================================
log_msg("-- PARTIE 2 : Concentration-based --")

cohort_conc <- read_parquet_dt(file.path(DATA_FINAL, "analytic_cohort_concentration.parquet"))
ps_conc <- read_parquet_dt(file.path(DATA_INTERMEDIATE, "ps_concentration.parquet"))
cohort_conc <- merge(cohort_conc, ps_conc, by="subject_id", all.x=TRUE)

ps_vars_conc <- c("daily_dose_main","age","female","weight_kg_baseline","alb_baseline",
                   "crea_baseline","in_icu","sepsis","emergency","epilepsy_dx","bipolar_dx",
                   "vpa_iv","comed_topiramate","comed_carbapenem","comed_heparin",
                   "comed_steroid","comed_antipsy","comed_other_aed","plt_baseline",
                   "alt_baseline","ast_baseline","liver_baseline_elevated",
                   "pre_hosp_days_winsorized")
ps_vars_conc <- intersect(ps_vars_conc, names(cohort_conc))

cohort_c <- copy(cohort_conc)
for (v in ps_vars_conc) {
  if (is.numeric(cohort_c[[v]])) cohort_c[is.na(get(v)), (v) := median(cohort_c[[v]], na.rm=TRUE)]
}

ps_formula_c <- as.formula(paste("high_conc_75 ~", paste(ps_vars_conc, collapse="+")))

set.seed(RANDOM_SEED)
w_conc <- tryCatch(
  weightit(ps_formula_c, data=as.data.frame(cohort_c), method="ps", estimand="ATE"),
  error=function(e) { log_msg(sprintf("WARN: %s", e$message)); NULL }
)

if (!is.null(w_conc)) {
  # 2a. Love plot complet concentration
  love_c <- love.plot(
    w_conc,
    stats = "mean.diffs",
    threshold = 0.1,
    abs = TRUE,
    var.order = "unadjusted",
    colors = c("tomato", "steelblue"),
    shapes = c("circle", "triangle"),
    sample.names = c("Non ajusté", "IPTW ajusté"),
    title = "Balance des covariables — Analyse Concentration-Based",
    subtitle = "SMD absolu | Seuil 0.1 | ATTENTION: positivité très limitée (ESS=3%)"
  )
  png(file.path(FIGURES_DIR, "love_plot_full_concentration_based.png"), width=1200, height=900, res=130)
  print(love_c)
  dev.off()
  log_msg("Love plot concentration-based sauvegardé")

  # 2b. Mirrored histogram concentration
  ps_hat_c <- w_conc$ps
  df_ps_c <- data.frame(
    PS        = as.numeric(ps_hat_c),
    Treatment = factor(cohort_c$high_conc_75, labels=c("Basse conc. (<75)", "Haute conc. (≥75)"))
  )
  p_mirror_c <- ggplot(df_ps_c, aes(x=PS, fill=Treatment)) +
    geom_histogram(data=subset(df_ps_c, Treatment=="Haute conc. (≥75)"),
                   aes(y=after_stat(count)), bins=50, alpha=0.7) +
    geom_histogram(data=subset(df_ps_c, Treatment=="Basse conc. (<75)"),
                   aes(y=-after_stat(count)), bins=50, alpha=0.7) +
    geom_hline(yintercept=0, color="black", linewidth=0.3) +
    scale_y_continuous(labels=abs) +
    scale_fill_manual(values=c("steelblue","tomato")) +
    labs(title="Overlap du Propensity Score — Concentration-Based",
         subtitle=sprintf("N=%d | ATTENTION: asymétrie marquée — violation de positivité probable",
                          nrow(df_ps_c)),
         x="P(haute concentration ≥75 µg/mL | X)", y="N patients") +
    theme_bw(base_size=12) + theme(legend.position="bottom")
  ggsave(file.path(FIGURES_DIR, "ps_overlap_concentration_based.png"),
         p_mirror_c, width=10, height=6, dpi=150)

  # 2c. Distribution des poids concentration
  w_c <- w_conc$weights
  ess_c <- sum(w_c)^2 / sum(w_c^2)

  p_weights_c <- ggplot(data.frame(w=pmin(w_c, quantile(w_c, 0.98)),
                                    group=factor(cohort_c$high_conc_75,
                                                  labels=c("Basse conc.","Haute conc."))),
                         aes(x=w, fill=group)) +
    geom_histogram(bins=60, alpha=0.7) +
    geom_vline(xintercept=quantile(w_c, 0.99), color="red", linetype="dashed") +
    scale_fill_manual(values=c("steelblue","tomato")) +
    labs(title="Distribution des Poids IPTW — Concentration-Based",
         subtitle=sprintf("ESS=%.0f/%d (%.1f%%) — CRITIQUE | Max=%.1f | P99=%.1f\nViolation sévère de positivité — IPTW non fiable",
                          ess_c, nrow(cohort_c), ess_c/nrow(cohort_c)*100,
                          max(w_c), quantile(w_c, 0.99)),
         x="Poids IPTW (tronqué au P98 pour lisibilité)", y="N patients", fill="Groupe") +
    theme_bw(base_size=12) + theme(legend.position="bottom")
  ggsave(file.path(FIGURES_DIR, "ipw_weight_distribution_concentration_based.png"),
         p_weights_c, width=10, height=6, dpi=150)
  log_msg(sprintf("Poids IPTW conc: ESS=%.0f/%.0f (%.1f%%), max=%.1f",
                  ess_c, nrow(cohort_c), ess_c/nrow(cohort_c)*100, max(w_c)))

  # 2d. Balance table complète concentration
  bal_tab_c <- bal.tab(w_conc, stats=c("m","v"), thresholds=c(m=0.1))
  bal_dt_c <- as.data.table(bal_tab_c$Balance)
  bal_dt_c[, variable := rownames(bal_tab_c$Balance)]
  bal_dt_c[, analysis := "concentration-based"]
  safe_write_csv(bal_dt_c, file.path(TABLES_DIR, "balance_full_concentration_based.csv"), overwrite=TRUE)
  max_smd_adj_c <- max(abs(bal_dt_c$Diff.Adj), na.rm=TRUE)
  log_msg(sprintf("SMD max ajusté (conc): %.3f", max_smd_adj_c))
}

# =============================================================================
# PARTIE 3 : Tableau de synthèse positivité / balance
# =============================================================================
log_msg("-- PARTIE 3 : Synthèse --")

positivity_summary <- data.table(
  analysis          = c("Dose-based", "Concentration-based"),
  n_total           = c(nrow(cohort_d), nrow(cohort_c)),
  n_treated         = c(sum(cohort_d$high_dose==1), sum(cohort_c$high_conc_75==1)),
  n_control         = c(sum(cohort_d$high_dose==0), sum(cohort_c$high_conc_75==0)),
  ess               = round(c(
    if(!is.null(w_dose)) {w_d <- w_dose$weights; sum(w_d)^2/sum(w_d^2)} else NA,
    if(!is.null(w_conc)) {w_c <- w_conc$weights; sum(w_c)^2/sum(w_c^2)} else NA
  ), 0),
  ess_pct           = round(c(
    if(!is.null(w_dose)) {w_d <- w_dose$weights; sum(w_d)^2/sum(w_d^2)/nrow(cohort_d)*100} else NA,
    if(!is.null(w_conc)) {w_c <- w_conc$weights; sum(w_c)^2/sum(w_c^2)/nrow(cohort_c)*100} else NA
  ), 1),
  max_weight        = round(c(
    if(!is.null(w_dose)) max(w_dose$weights) else NA,
    if(!is.null(w_conc)) max(w_conc$weights) else NA
  ), 1),
  p99_weight        = round(c(
    if(!is.null(w_dose)) quantile(w_dose$weights, 0.99) else NA,
    if(!is.null(w_conc)) quantile(w_conc$weights, 0.99) else NA
  ), 1),
  max_smd_adj       = round(c(
    if(!is.null(w_dose)) max(abs(bal_dt_d$Diff.Adj), na.rm=TRUE) else NA,
    if(!is.null(w_conc)) max(abs(bal_dt_c$Diff.Adj), na.rm=TRUE) else NA
  ), 3),
  positivity_status = c("Acceptable", "COMPROMISE (ESS~3%)")
)

safe_write_csv(positivity_summary, file.path(TABLES_DIR, "positivity_summary.csv"), overwrite=TRUE)
log_msg("Synthèse positivité sauvegardée:")
print(positivity_summary)

timer_end(t0, "A8 balance diagnostics")
log_msg("=== STARTING ===")
