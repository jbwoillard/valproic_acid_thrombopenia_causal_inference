# =============================================================================
# describe_population.R
# Step 3 — Descriptive statistics, Table 1, outcome summaries
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(ggplot2)
library(patchwork)
library(gtsummary)
library(cobalt)
library(arrow)

log_msg("=== ÉTAPE 7 : ANALYSE DESCRIPTIVE ===")
t0 <- timer_start()

cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
log_msg(sprintf("Cohorte: %d patients", nrow(cohort)))


# =============================================================================
# 1. Table 1 — Avant ajustement seul et split par groupe de dose
# =============================================================================
log_msg("1. Table 1...")

library(data.table)
library(dplyr)
library(gtsummary)
library(gt)

table1_vars <- c(
  "age", "female", "race_simple", "emergency",
  "in_icu", "sepsis", "epilepsy_dx", "bipolar_dx",
  "plt_baseline", "alb_baseline", "crea_baseline",
  "alt_baseline", "ast_baseline",
  "liver_baseline_elevated", "alb_baseline_low", "crea_aki",
  "comed_topiramate", "comed_carbapenem", "comed_heparin",
  "comed_steroid", "comed_other_aed", "vpa_iv",
  "daily_dose_main", "severity_score", "pre_hosp_days_winsorized"
)

table1_vars_ok <- table1_vars[table1_vars %in% names(cohort)]

# Copie de travail
cohort_tbl <- as.data.frame(cohort[, c(table1_vars_ok, "high_dose"), with = FALSE])

# Variable groupe lisible
cohort_tbl$high_dose <- factor(
  cohort_tbl$high_dose,
  levels = c(0, 1),
  labels = c("Low dose (≤1000 mg/day)", "High dose (>1000 mg/day)")
)

# Forcer les variables binaires en facteur pour affichage n (%)
binary_vars <- c(
  "female", "emergency", "in_icu", "sepsis", "epilepsy_dx", "bipolar_dx",
  "liver_baseline_elevated", "alb_baseline_low", "crea_aki",
  "comed_topiramate", "comed_carbapenem", "comed_heparin",
  "comed_steroid", "comed_other_aed", "vpa_iv"
)
binary_vars_ok <- intersect(binary_vars, names(cohort_tbl))

cohort_tbl[binary_vars_ok] <- lapply(cohort_tbl[binary_vars_ok], function(x) {
  factor(x, levels = c(0, 1), labels = c("No", "Yes"))
})

# Table gtsummary
tbl1 <- tryCatch({
  tbl_summary(
    cohort_tbl,
    by = high_dose,
    missing = "ifany",
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age (years)",
      female ~ "Female sex",
      race_simple ~ "Race/ethnicity",
      emergency ~ "Emergency admission",
      in_icu ~ "ICU stay",
      sepsis ~ "Sepsis (ICD)",
      epilepsy_dx ~ "Epilepsy (ICD)",
      bipolar_dx ~ "Bipolar disorder (ICD)",
      plt_baseline ~ "Baseline platelets (×10^9/L)",
      alb_baseline ~ "Baseline albumin (g/dL)",
      crea_baseline ~ "Baseline creatinine (mg/dL)",
      alt_baseline ~ "Baseline ALT (U/L)",
      ast_baseline ~ "Baseline AST (U/L)",
      liver_baseline_elevated ~ "Pre-existing liver abnormality",
      alb_baseline_low ~ "Hypoalbuminaemia (<3 g/dL)",
      crea_aki ~ "Baseline renal impairment",
      comed_topiramate ~ "Concomitant topiramate",
      comed_carbapenem ~ "Concomitant carbapenem",
      comed_heparin ~ "Concomitant heparin",
      comed_steroid ~ "Concomitant steroid",
      comed_other_aed ~ "Other antiseizure drug",
      vpa_iv ~ "Intravenous valproate",
      daily_dose_main ~ "Valproate daily dose (mg/day)",
      severity_score ~ "Approximate severity score",
      pre_hosp_days_winsorized ~ "Days from admission to valproate"
    )
  ) |>
    add_overall(last = FALSE) |>
    add_p() |>
    add_n() |>
    bold_labels()
}, error = function(e) {
  log_msg(sprintf("Table 1 simplifiée (erreur gtsummary: %s)", e$message), level = "WARN")
  NULL
})

if (!is.null(tbl1)) {
  # HTML
  tryCatch({
    gt_tbl <- as_gt(tbl1)
    gt::gtsave(gt_tbl, file.path(TABLES_DIR, "table1_unadjusted.html"))
    log_msg("Table 1 HTML sauvegardée")
  }, error = function(e) log_msg(sprintf("Impossible de sauvegarder en HTML: %s", e$message), level = "WARN"))
  
  # CSV exploitable pour Word / Excel
  tryCatch({
    t1_body <- as_tibble(tbl1$table_body)
    
    # garder uniquement les colonnes utiles
    keep_cols <- intersect(
      c("label", "stat_0", "stat_1", "stat_2", "p.value"),
      names(t1_body)
    )
    
    t1_csv <- t1_body[, keep_cols]
    
    # renommer proprement
    names(t1_csv) <- c(
      "Characteristic",
      "Overall",
      "Low dose (≤1000 mg/day)",
      "High dose (>1000 mg/day)",
      "P_value"
    )[seq_along(keep_cols)]
    
    safe_write_csv(as.data.table(t1_csv), file.path(TABLES_DIR, "table1_unadjusted.csv"), overwrite = TRUE)
    log_msg("Table 1 CSV par groupe sauvegardée")
  }, error = function(e) {
    log_msg(sprintf("Impossible de sauvegarder Table 1 CSV: %s", e$message), level = "WARN")
  })
}

# =============================================================================
# SMD non ajustés
# =============================================================================
# Fonction SMD simple
calc_smd <- function(x, g) {
  # g doit être binaire 0/1
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  
  if (length(unique(g)) != 2) return(NA_real_)
  
  # Numérique
  if (is.numeric(x)) {
    m1 <- mean(x[g == 1], na.rm = TRUE)
    m0 <- mean(x[g == 0], na.rm = TRUE)
    s1 <- stats::var(x[g == 1], na.rm = TRUE)
    s0 <- stats::var(x[g == 0], na.rm = TRUE)
    sp <- sqrt((s1 + s0) / 2)
    return(ifelse(is.finite(sp) && sp > 0, abs((m1 - m0) / sp), NA_real_))
  }
  
  # Binaire / facteur à 2 niveaux
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    x2 <- as.factor(x)
    if (nlevels(x2) == 2) {
      p1 <- mean(x2[g == 1] == levels(x2)[2], na.rm = TRUE)
      p0 <- mean(x2[g == 0] == levels(x2)[2], na.rm = TRUE)
      p  <- (p1 + p0) / 2
      denom <- sqrt(p * (1 - p))
      return(ifelse(is.finite(denom) && denom > 0, abs((p1 - p0) / denom), NA_real_))
    }
  }
  
  return(NA_real_)
}

smd_vars <- intersect(table1_vars_ok, names(cohort))

smd_results <- sapply(smd_vars, function(v) {
  calc_smd(cohort[[v]], cohort$high_dose)
})

smd_dt <- data.table(
  variable = names(smd_results),
  smd_unadjusted = round(as.numeric(smd_results), 3)
)

safe_write_csv(smd_dt, file.path(TABLES_DIR, "smd_unadjusted.csv"), overwrite = TRUE)
log_msg("SMD calculés:")
print(smd_dt)

# =============================================================================
# 2. Distribution de l'exposition
# =============================================================================
log_msg("2. Distribution de l'exposition...")

p_dose_dist <- ggplot(cohort, aes(x=daily_dose_main, fill=factor(high_dose))) +
  geom_histogram(binwidth=100, alpha=0.7, position="identity") +
  geom_vline(xintercept=1000, color="red", linetype="dashed", linewidth=1) +
  scale_fill_manual(values=c("steelblue","tomato"),
                    labels=c("Faible dose (≤1000)","Haute dose (>1000)")) +
  labs(title="Distribution de la dose journalière de valproate",
       subtitle="Ligne rouge = seuil 1000 mg/j (définition principale de l'exposition)",
       x="Dose journalière (mg/j)", y="N patients", fill="Groupe") +
  theme_bw(base_size=12)

p_dose_box <- ggplot(cohort, aes(x=factor(high_dose), y=daily_dose_main,
                                   fill=factor(high_dose))) +
  geom_boxplot(alpha=0.7, outlier.size=0.5) +
  scale_fill_manual(values=c("steelblue","tomato"),
                    labels=c("Faible","Haute")) +
  labs(x="Groupe de dose", y="Dose journalière (mg/j)", fill="Groupe") +
  theme_bw(base_size=12) + theme(legend.position="none")

p_dose_combined <- p_dose_dist / p_dose_box + plot_layout(heights=c(3,1))
ggsave(file.path(FIGURES_DIR, "descriptive_dose_distribution.png"),
       p_dose_combined, width=10, height=8, dpi=150)

# =============================================================================
# 3. Missingness heatmap
# =============================================================================
log_msg("3. Heatmap de données manquantes...")

miss_vars <- c("plt_baseline","alb_baseline","crea_baseline","alt_baseline",
               "ast_baseline","bili_baseline","nh3_baseline","weight_kg_baseline",
               "vpa_level_mean_wk1")
miss_vars_ok <- intersect(miss_vars, names(cohort))

miss_long <- melt(
  cohort[, c("subject_id", miss_vars_ok), with=FALSE],
  id.vars="subject_id",
  variable.name="variable",
  value.name="value"
)
miss_long[, missing := as.integer(is.na(value))]

miss_by_var <- miss_long[, .(pct_missing = mean(missing)*100), by=variable]

p_miss <- ggplot(miss_by_var, aes(x=reorder(variable, pct_missing), y=pct_missing)) +
  geom_col(fill="steelblue", alpha=0.8) +
  geom_hline(yintercept=40, color="orange", linetype="dashed") +
  geom_hline(yintercept=20, color="green", linetype="dashed") +
  coord_flip() +
  labs(title="Données manquantes par variable",
       subtitle="Lignes: rouge=40%, vert=20%",
       x="Variable", y="% manquant") +
  theme_bw(base_size=12)

ggsave(file.path(FIGURES_DIR, "descriptive_missingness.png"),
       p_miss, width=10, height=6, dpi=150)

# =============================================================================
# 4. Distribution des outcomes
# =============================================================================
log_msg("4. Distribution des outcomes...")

outcome_summary <- cohort[, .(
  n = .N,
  thrombo_100  = sum(outcome_thrombo_100, na.rm=T),
  pct_thrombo  = mean(outcome_thrombo_100, na.rm=T)*100,
  thrombo_50   = sum(outcome_thrombo_50, na.rm=T),
  pct_thrombo50 = mean(outcome_thrombo_50, na.rm=T)*100,
  thrombo_rel30 = sum(outcome_thrombo_rel30, na.rm=T),
  pct_rel30    = mean(outcome_thrombo_rel30, na.rm=T)*100,
  hepato       = sum(outcome_hepato_3xULN, na.rm=T),
  pct_hepato   = mean(outcome_hepato_3xULN, na.rm=T)*100,
  hypernh3     = sum(outcome_hypernh3, na.rm=T),
  pct_nh3      = mean(outcome_hypernh3, na.rm=T)*100,
  composite    = sum(outcome_composite, na.rm=T),
  pct_comp     = mean(outcome_composite, na.rm=T)*100
), by=high_dose]
print(outcome_summary)
safe_write_csv(outcome_summary, file.path(TABLES_DIR, "outcome_summary_by_group.csv"), overwrite=TRUE)

# Figure crude outcome rates
outcome_long <- melt(
  outcome_summary,
  id.vars="high_dose",
  measure.vars=c("pct_thrombo","pct_thrombo50","pct_rel30","pct_hepato","pct_nh3"),
  variable.name="outcome", value.name="pct"
)
outcome_long[, outcome_label := fcase(
  outcome == "pct_thrombo",   "Thrombopénie\n< 100 G/L",
  outcome == "pct_thrombo50", "Thrombopénie\n< 50 G/L",
  outcome == "pct_rel30",     "Baisse plt\n> 30%",
  outcome == "pct_hepato",    "Hépatotoxicité\nALT/AST>3xULN",
  outcome == "pct_nh3",       "Hyperammoniémie\nNH3>55"
)]
outcome_long[, group_label := factor(high_dose, labels=c("Faible dose","Haute dose"))]

p_outcomes <- ggplot(outcome_long, aes(x=outcome_label, y=pct, fill=group_label)) +
  geom_col(position="dodge", alpha=0.8) +
  scale_fill_manual(values=c("steelblue","tomato")) +
  labs(title="Taux d'outcomes par groupe de dose (non ajusté)",
       x="Outcome", y="Prévalence (%)", fill="Groupe") +
  theme_bw(base_size=12)

ggsave(file.path(FIGURES_DIR, "descriptive_outcomes_crude.png"),
       p_outcomes, width=12, height=6, dpi=150)

# =============================================================================
# 5. Distribution des plaquettes baseline vs plaquettes minimales en suivi
# =============================================================================
log_msg("5. Distribution des plaquettes...")

plt_data <- cohort[!is.na(plt_baseline) & !is.na(plt_min_fu), .(
  plt_baseline, plt_min_fu, high_dose,
  delta_plt = plt_baseline - plt_min_fu
)]

p_plt_scatter <- ggplot(plt_data, aes(x=plt_baseline, y=plt_min_fu, color=factor(high_dose))) +
  geom_point(alpha=0.3, size=0.8) +
  geom_abline(slope=1, intercept=0, color="black", linetype="dashed") +
  geom_hline(yintercept=100, color="red", linetype="dotted") +
  scale_color_manual(values=c("steelblue","tomato"), labels=c("Faible dose","Haute dose")) +
  labs(title="Plaquettes baseline vs minimum sur 30j de suivi",
       subtitle="Ligne rouge: seuil thrombopénie; ligne noire: identité",
       x="Plaquettes baseline (G/L)", y="Plaquettes min suivi (G/L)", color="Groupe") +
  xlim(0, 600) + ylim(0, 600) +
  theme_bw(base_size=12)

ggsave(file.path(FIGURES_DIR, "descriptive_platelets_scatter.png"),
       p_plt_scatter, width=8, height=7, dpi=150)

# =============================================================================
# 6. Diagnostics de positivity
# =============================================================================
log_msg("6. Diagnostics de positivity préliminaires...")

# Croisement ICU x haute dose
tbl_icu_dose <- table(cohort$in_icu, cohort$high_dose,
                       dnn=c("ICU","Haute dose"))
log_msg("Croisement ICU x Haute dose:")
print(tbl_icu_dose)
cat("Proportions par ICU:\n")
print(prop.table(tbl_icu_dose, 1))

# Vérification overlap : sévérité vs dose
p_overlap <- ggplot(cohort, aes(x=severity_score, fill=factor(high_dose))) +
  geom_histogram(binwidth=1, position="identity", alpha=0.6) +
  scale_fill_manual(values=c("steelblue","tomato"), labels=c("Faible dose","Haute dose")) +
  labs(title="Overlap : Score de sévérité par groupe de dose",
       subtitle="Overlap nécessaire pour les méthodes causales",
       x="Score de sévérité approché", y="N patients", fill="Groupe") +
  theme_bw(base_size=12)

ggsave(file.path(FIGURES_DIR, "descriptive_overlap_severity.png"),
       p_overlap, width=9, height=6, dpi=150)

timer_end(t0, "Étape 7 totale")
log_msg("Figures sauvegardées dans:", FIGURES_DIR)

# Rapport minimal
writeLines(sprintf(
"# Analyse Descriptive — Étude VPA/MIMIC-IV
**Date :** %s

## Population globale
- N = 2,640 patients
- Haute dose (>1000 mg/j) : 884 (33.5%%)
- Faible dose (<=1000 mg/j) : 1,756 (66.5%%)

## Taux d'outcomes par groupe
%s

## Données manquantes critiques
- NH3 baseline : 97.5%% manquant (mesure peu disponible)
- Albumine baseline : 51.6%% manquant
- Poids baseline : 59.9%% manquant
- ALT/AST baseline : 37-38%% manquant

## Figures produites
- `descriptive_dose_distribution.png`
- `descriptive_missingness.png`
- `descriptive_outcomes_crude.png`
- `descriptive_platelets_scatter.png`
- `descriptive_overlap_severity.png`
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(apply(outcome_summary, 1, function(r)
    sprintf("- %s: thrombo=%.1f%%, hepato=%.1f%%",
            ifelse(as.integer(r["high_dose"])==1, "Haute dose", "Faible dose"),
            as.numeric(r["pct_thrombo"]),
            as.numeric(r["pct_hepato"]))),
    collapse="\n")
), file.path(REPORTS_DIR, "07_descriptive_analysis.md"))

log_msg("=== ÉTAPE 7 TERMINÉE ===")
