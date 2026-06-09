# =============================================================================
# build_concentration_subcohort.R
# Conc-2 — Concentration-based subcohort construction
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(arrow)
library(ggplot2)
library(patchwork)

log_msg("=== ÉTAPE 8B.4-8B.5 : COHORTE ET DESCRIPTIF CONCENTRATION-BASED ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR, DIAGNOSTICS_DIR)

# =============================================================================
# 8B.4 — Construction de la sous-cohorte monitorée
# =============================================================================
log_msg("8B.4 — Construction de la sous-cohorte monitorée...")

cohort_full <- read_parquet_dt(file.path(DATA_INTERMEDIATE, "cohort_with_monitoring_flag.parquet"))
vpa_early   <- read_parquet_dt(file.path(DATA_INTERMEDIATE, "vpa_levels_wk1.parquet"))

log_msg(sprintf("Cohorte complète: %d patients", nrow(cohort_full)))
log_msg(sprintf("VPA levels J0-J7: %d mesures, %d patients",
                nrow(vpa_early), uniqueN(vpa_early$subject_id)))

# Exclure valeurs aberrantes > 500 µg/mL
n_outlier <- sum(vpa_early$value_ugml > 500)
vpa_early  <- vpa_early[value_ugml <= 500]
log_msg(sprintf("Exclusions > 500 µg/mL: %d mesures exclues", n_outlier))

# Calcul de l'exposition concentration-based par patient
conc_by_pt <- vpa_early[, .(
  first_conc   = value_ugml[which.min(days_from_tz)],
  first_conc_t = days_from_tz[which.min(days_from_tz)],
  max_conc     = max(value_ugml),
  mean_conc    = mean(value_ugml),
  last_conc    = value_ugml[which.max(days_from_tz)],
  n_levels     = .N
), by=subject_id]

log_msg(sprintf("Patients avec concentration calculée: %d", nrow(conc_by_pt)))

# Définitions de l'exposition binaire
conc_by_pt[, high_conc_75  := as.integer(first_conc >= 75)]   # Seuil principal
conc_by_pt[, high_conc_50  := as.integer(first_conc >= 50)]   # Sensibilité 1
conc_by_pt[, high_conc_100 := as.integer(first_conc >= 100)]  # Sensibilité 2

# Quartiles de concentration
q_conc <- quantile(conc_by_pt$first_conc, c(0.25, 0.5, 0.75), na.rm=TRUE)
conc_by_pt[, conc_quartile := cut(first_conc,
  breaks = c(-Inf, q_conc[1], q_conc[2], q_conc[3], Inf),
  labels = c("Q1 (bas)", "Q2", "Q3", "Q4 (haut)"),
  include.lowest = TRUE
)]

log_msg(sprintf("Seuils quartiles concentration: Q1=%.1f, Q2=%.1f, Q3=%.1f µg/mL",
                q_conc[1], q_conc[2], q_conc[3]))
log_msg(sprintf("Haute concentration (≥75): %d (%.1f%%)",
                sum(conc_by_pt$high_conc_75), mean(conc_by_pt$high_conc_75)*100))

# Fusion avec la cohorte principale
cohort_conc <- merge(
  cohort_full[monitored_early == 1],
  conc_by_pt,
  by="subject_id",
  all.x=FALSE
)

log_msg(sprintf("Sous-cohorte monitorée après fusion: %d patients", nrow(cohort_conc)))

# Vérification des outcomes disponibles
for (oc in c("outcome_thrombo_100","outcome_thrombo_50","outcome_hepato_3xULN",
             "outcome_hypernh3","outcome_composite")) {
  n_ev <- sum(cohort_conc[[oc]], na.rm=TRUE)
  log_msg(sprintf("  %s: %d événements (%.1f%%)", oc, n_ev, n_ev/nrow(cohort_conc)*100))
}

# Sauvegarde
arrow::write_parquet(as.data.frame(cohort_conc),
                      file.path(DATA_FINAL, "analytic_cohort_concentration.parquet"))
log_msg(sprintf("Cohorte sauvegardée: data_final/analytic_cohort_concentration.parquet"))

# Cohort flow
cohort_flow <- data.table(
  step = c(
    "1. Cohorte dose-based principale",
    "2. Patients avec dosage VPA dans J0-J7 (monitored_early=1)",
    "3. Exclusion concentrations > 500 µg/mL",
    "4. Fusion avec concentrations calculées"
  ),
  n_patients = c(
    nrow(cohort_full),
    sum(cohort_full$monitored_early),
    sum(cohort_full$monitored_early),
    nrow(cohort_conc)
  ),
  exclusions = c(NA, nrow(cohort_full) - sum(cohort_full$monitored_early),
                 n_outlier, 0)
)
safe_write_csv(cohort_flow, file.path(TABLES_DIR, "cohort_flow_concentration.csv"), overwrite=TRUE)

# =============================================================================
# 8B.5 — Analyse descriptive de la sous-cohorte
# =============================================================================
log_msg("8B.5 — Analyse descriptive...")

library(data.table)

# -----------------------------------------------------------------------------
# Variables à décrire
# -----------------------------------------------------------------------------
desc_vars <- c(
  "age", "female", "in_icu", "sepsis", "emergency", "epilepsy_dx", "bipolar_dx",
  "vpa_iv", "daily_dose_main", "plt_baseline", "alb_baseline", "crea_baseline",
  "alt_baseline", "ast_baseline", "n_levels", "first_conc_t",
  "outcome_thrombo_100", "outcome_thrombo_50", "outcome_hepato_3xULN",
  "outcome_hypernh3", "outcome_composite"
)
desc_vars <- intersect(desc_vars, names(cohort_conc))

# Variables binaires explicites
binary_vars <- intersect(c(
  "female", "in_icu", "sepsis", "emergency", "epilepsy_dx", "bipolar_dx",
  "vpa_iv", "outcome_thrombo_100", "outcome_thrombo_50",
  "outcome_hepato_3xULN", "outcome_hypernh3", "outcome_composite"
), desc_vars)

# Variables continues explicites
continuous_vars <- setdiff(desc_vars, binary_vars)

# -----------------------------------------------------------------------------
# Libellés propres
# -----------------------------------------------------------------------------
var_labels <- c(
  age = "Age, years",
  female = "Female sex",
  in_icu = "ICU stay",
  sepsis = "Sepsis (ICD)",
  emergency = "Emergency admission",
  epilepsy_dx = "Epilepsy (ICD)",
  bipolar_dx = "Bipolar disorder (ICD)",
  vpa_iv = "Intravenous valproate",
  daily_dose_main = "Daily valproate dose, mg/day",
  plt_baseline = "Baseline platelets, ×10^9/L",
  alb_baseline = "Baseline albumin, g/dL",
  crea_baseline = "Baseline creatinine, mg/dL",
  alt_baseline = "Baseline ALT, U/L",
  ast_baseline = "Baseline AST, U/L",
  n_levels = "Number of concentration measurements",
  first_conc_t = "Time to first concentration, days",
  outcome_thrombo_100 = "Incident thrombocytopenia <100 ×10^9/L",
  outcome_thrombo_50 = "Severe thrombocytopenia <50 ×10^9/L",
  outcome_hepato_3xULN = "Hepatotoxicity >3×ULN",
  outcome_hypernh3 = "Hyperammonaemia >55 µmol/L",
  outcome_composite = "Composite outcome"
)

# -----------------------------------------------------------------------------
# Noms lisibles des groupes
# -----------------------------------------------------------------------------
cohort_conc_copy <- copy(cohort_conc)
cohort_conc_copy[, conc_group := fifelse(
  high_conc_75 == 1,
  "Higher concentration (≥75 µg/mL)",
  "Lower concentration (<75 µg/mL)"
)]

# -----------------------------------------------------------------------------
# Fonctions de résumé
# -----------------------------------------------------------------------------
fmt_cont <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  sprintf(
    "%.1f (%.1f–%.1f)",
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

fmt_bin <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  n1 <- sum(x == 1, na.rm = TRUE)
  denom <- sum(!is.na(x))
  sprintf("%d (%.1f%%)", n1, 100 * n1 / denom)
}

# -----------------------------------------------------------------------------
# Construction de la table S1
# -----------------------------------------------------------------------------
rows_list <- list()

for (v in desc_vars) {
  overall <- if (v %in% continuous_vars) {
    fmt_cont(cohort_conc_copy[[v]])
  } else {
    fmt_bin(cohort_conc_copy[[v]])
  }
  
  low <- if (v %in% continuous_vars) {
    fmt_cont(cohort_conc_copy[conc_group == "Lower concentration (<75 µg/mL)"][[v]])
  } else {
    fmt_bin(cohort_conc_copy[conc_group == "Lower concentration (<75 µg/mL)"][[v]])
  }
  
  high <- if (v %in% continuous_vars) {
    fmt_cont(cohort_conc_copy[conc_group == "Higher concentration (≥75 µg/mL)"][[v]])
  } else {
    fmt_bin(cohort_conc_copy[conc_group == "Higher concentration (≥75 µg/mL)"][[v]])
  }
  
  rows_list[[v]] <- data.table(
    Characteristic = unname(var_labels[v]),
    Overall = overall,
    `Lower concentration (<75 µg/mL)` = low,
    `Higher concentration (≥75 µg/mL)` = high
  )
}

table1_conc <- rbindlist(rows_list, use.names = TRUE, fill = TRUE)

safe_write_csv(
  table1_conc,
  file.path(TABLES_DIR, "table1_concentration_unadjusted.csv"),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# SMD par groupe de concentration
# -----------------------------------------------------------------------------
calc_smd <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  
  if (length(unique(g)) != 2) return(NA_real_)
  
  if (is.numeric(x) && !all(x %in% c(0, 1), na.rm = TRUE)) {
    m1 <- mean(x[g == 1], na.rm = TRUE)
    m0 <- mean(x[g == 0], na.rm = TRUE)
    s1 <- stats::var(x[g == 1], na.rm = TRUE)
    s0 <- stats::var(x[g == 0], na.rm = TRUE)
    sp <- sqrt((s1 + s0) / 2)
    return(ifelse(is.finite(sp) && sp > 0, abs((m1 - m0) / sp), NA_real_))
  }
  
  # Variables binaires 0/1
  p1 <- mean(x[g == 1] == 1, na.rm = TRUE)
  p0 <- mean(x[g == 0] == 1, na.rm = TRUE)
  p <- (p1 + p0) / 2
  denom <- sqrt(p * (1 - p))
  return(ifelse(is.finite(denom) && denom > 0, abs((p1 - p0) / denom), NA_real_))
}

smd_vars_conc <- intersect(
  c("age", "female", "in_icu", "sepsis", "emergency", "epilepsy_dx", "bipolar_dx",
    "vpa_iv", "daily_dose_main", "plt_baseline", "alb_baseline", "crea_baseline",
    "alt_baseline", "ast_baseline", "n_levels", "first_conc_t"),
  names(cohort_conc_copy)
)

smd_conc <- sapply(smd_vars_conc, function(v) {
  tryCatch(calc_smd(cohort_conc_copy[[v]], cohort_conc_copy$high_conc_75),
           error = function(e) NA_real_)
})

smd_conc_dt <- data.table(
  variable = names(smd_conc),
  label = unname(var_labels[names(smd_conc)]),
  smd = round(as.numeric(smd_conc), 3)
)

safe_write_csv(
  smd_conc_dt,
  file.path(TABLES_DIR, "smd_concentration_groups.csv"),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# SMD monitored vs non-monitored
# -----------------------------------------------------------------------------
cohort_full[, monitored_flag := monitored_early]

smd_vars_monit <- intersect(
  c("age", "female", "in_icu", "sepsis", "epilepsy_dx", "vpa_iv", "daily_dose_main",
    "plt_baseline", "alb_baseline", "crea_baseline", "high_dose"),
  names(cohort_full)
)

smd_monit <- sapply(smd_vars_monit, function(v) {
  tryCatch(calc_smd(cohort_full[[v]], cohort_full$monitored_flag),
           error = function(e) NA_real_)
})

smd_monit_dt <- data.table(
  variable = names(smd_monit),
  smd = round(as.numeric(smd_monit), 3)
)

safe_write_csv(
  smd_monit_dt,
  file.path(TABLES_DIR, "smd_monitored_vs_nonmonitored.csv"),
  overwrite = TRUE
)

log_msg("Table S1 et SMD concentration groups générés")


# =============================================================================
# Figure 5 — Exploratory concentration-based analysis
# Publication-ready version for BJCP
# =============================================================================

library(ggplot2)
library(patchwork)
library(data.table)

# Optional clean theme
theme_bjcp <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

# ------------------------------------------------------------------
# Panel A — Distribution of first valproate concentration by dose group
# ------------------------------------------------------------------
p5a <- ggplot(
  cohort_conc[first_conc <= 300],
  aes(x = first_conc, fill = factor(high_dose))
) +
  geom_histogram(
    bins = 35,
    alpha = 0.65,
    position = "identity",
    colour = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 75,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ factor(high_dose,
             levels = c(0, 1),
             labels = c("Lower dose (≤1000 mg/day)",
                        "Higher dose (>1000 mg/day)"))
  ) +
  scale_fill_manual(
    values = c("grey70", "grey30"),
    labels = c("Lower dose", "Higher dose")
  ) +
  labs(
    title = "A. Distribution of first valproate concentration",
    x = "First valproate concentration (µg/mL)",
    y = "Number of patients",
    fill = "Dose group"
  ) +
  coord_cartesian(xlim = c(0, 300)) +
  theme_bjcp() +
  theme(legend.position = "none")

# ------------------------------------------------------------------
# Panel B — Prescribed daily dose versus first concentration
# ------------------------------------------------------------------
rho <- suppressWarnings(
  cor(
    cohort_conc$daily_dose_main,
    cohort_conc$first_conc,
    method = "spearman",
    use = "complete.obs"
  )
)

p5b <- ggplot(
  cohort_conc[first_conc <= 300],
  aes(x = daily_dose_main, y = first_conc)
) +
  geom_point(
    alpha = 0.35,
    size = 1.3,
    colour = "grey30"
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = 75,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = 1000,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Spearman's \u03c1 = ", round(rho, 2)),
    hjust = 1.1, vjust = 1.5,
    size = 3.4
  ) +
  labs(
    title = "B. Prescribed dose versus first valproate concentration",
    x = "Prescribed daily dose (mg/day)",
    y = "First valproate concentration (µg/mL)"
  ) +
  coord_cartesian(xlim = c(0, max(cohort_conc$daily_dose_main, na.rm = TRUE)),
                  ylim = c(0, 300)) +
  theme_bjcp()

# ------------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------------
fig5 <- p5a / p5b + plot_layout(heights = c(1, 1))

ggsave(
  filename = file.path(FIGURES_DIR, "Figure_5_concentration_timing_reappraisal.png"),
  plot = fig5,
  width = 9,
  height = 8,
  dpi = 300
)

# Rapport 8B.4-8B.5
# =============================================================================
n_conc <- nrow(cohort_conc)
n_high <- sum(cohort_conc$high_conc_75)
n_low  <- n_conc - n_high

report_8b45 <- sprintf(
'# 8B Cohorte et Descriptif Concentration-Based
**Date :** %s

---

## 8B.4 — Construction de la sous-cohorte monitorée

### Flow de la cohorte

| Étape | N patients | Exclusions |
|-------|-----------|------------|
| Cohorte dose-based principale | %d | — |
| Patients avec dosage VPA J0-J7 | %d | %d exclus (non monitorés) |
| Exclusion concentrations > 500 µg/mL | %d | %d |
| **Sous-cohorte finale** | **%d** | — |

### Définition de l\'exposition

- **Haute concentration** : première mesure J0-J7 ≥ 75 µg/mL → **%d patients (%.1f%%)**
- **Basse concentration** : première mesure J0-J7 < 75 µg/mL → **%d patients (%.1f%%)**

Définitions alternatives disponibles :
- ≥ 50 µg/mL : %d (%.1f%%)
- ≥ 100 µg/mL : %d (%.1f%%)

### Disponibilité des outcomes

| Outcome | N événements | %% |
|---------|-------------|-----|
| Thrombopénie < 100 G/L | %d | %.1f%% |
| Thrombopénie < 50 G/L | %d | %.1f%% |
| Hépatotoxicité > 3xULN | %d | %.1f%% |
| Hyperammoniémie > 55 µmol/L | %d | %.1f%% |
| Outcome composite | %d | %.1f%% |

---

## 8B.5 — Analyse descriptive

### Distribution des concentrations

| Statistique | Valeur (µg/mL) |
|-------------|----------------|
| P10 | %.1f |
| P25 | %.1f |
| P50 (médiane) | %.1f |
| P75 | %.1f |
| P90 | %.1f |
| Proportion < 50 µg/mL | %.1f%% |
| Proportion 50-75 µg/mL | %.1f%% |
| Proportion 75-100 µg/mL | %.1f%% |
| Proportion > 100 µg/mL | %.1f%% |

### Délai initiation → premier dosage (médiane, IQR)

Médiane = %.1f jours [IQR: %.1f – %.1f jours]

### Corrélation dose-concentration

Corrélation de Spearman dose journalière vs première concentration : r = %.2f

### Biais de sélection confirmé

Les patients monitorés diffèrent des non-monitorés sur plusieurs caractéristiques clés
(voir vpa_tdm_monitoring_smd.csv). Ce biais de sélection est majeur et limite la
généralisabilité des résultats de la sous-cohorte.

### Limites

- Timing du prélèvement (trough/peak) inconnu — variabilité pharmacocinétique non contrôlée
- Sélection sur gravité, suspicion de toxicité, doses élevées
- Sous-cohorte non représentative de l\'ensemble des initiateurs incidents de VPA
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(cohort_full),
  sum(cohort_full$monitored_early), nrow(cohort_full) - sum(cohort_full$monitored_early),
  sum(cohort_full$monitored_early), n_outlier,
  n_conc,
  n_high, n_high/n_conc*100,
  n_low, n_low/n_conc*100,
  sum(cohort_conc$high_conc_50), mean(cohort_conc$high_conc_50)*100,
  sum(cohort_conc$high_conc_100), mean(cohort_conc$high_conc_100)*100,
  sum(cohort_conc$outcome_thrombo_100, na.rm=TRUE),
  mean(cohort_conc$outcome_thrombo_100, na.rm=TRUE)*100,
  sum(cohort_conc$outcome_thrombo_50, na.rm=TRUE),
  mean(cohort_conc$outcome_thrombo_50, na.rm=TRUE)*100,
  sum(cohort_conc$outcome_hepato_3xULN, na.rm=TRUE),
  mean(cohort_conc$outcome_hepato_3xULN, na.rm=TRUE)*100,
  sum(cohort_conc$outcome_hypernh3, na.rm=TRUE),
  mean(cohort_conc$outcome_hypernh3, na.rm=TRUE)*100,
  sum(cohort_conc$outcome_composite, na.rm=TRUE),
  mean(cohort_conc$outcome_composite, na.rm=TRUE)*100,
  quantile(cohort_conc$first_conc, 0.10, na.rm=TRUE),
  quantile(cohort_conc$first_conc, 0.25, na.rm=TRUE),
  median(cohort_conc$first_conc, na.rm=TRUE),
  quantile(cohort_conc$first_conc, 0.75, na.rm=TRUE),
  quantile(cohort_conc$first_conc, 0.90, na.rm=TRUE),
  mean(cohort_conc$first_conc < 50, na.rm=TRUE)*100,
  mean(cohort_conc$first_conc >= 50 & cohort_conc$first_conc < 75, na.rm=TRUE)*100,
  mean(cohort_conc$first_conc >= 75 & cohort_conc$first_conc < 100, na.rm=TRUE)*100,
  mean(cohort_conc$first_conc >= 100, na.rm=TRUE)*100,
  median(cohort_conc$first_conc_t, na.rm=TRUE),
  quantile(cohort_conc$first_conc_t, 0.25, na.rm=TRUE),
  quantile(cohort_conc$first_conc_t, 0.75, na.rm=TRUE),
  cor(cohort_conc$daily_dose_main, cohort_conc$first_conc,
      method="spearman", use="complete.obs")
)

writeLines(report_8b45, file.path(REPORTS_DIR, "08B_cohort_concentration.md"))
log_msg("Rapport 8B.4-8B.5 écrit")

timer_end(t0, "8B.4-8B.5 cohorte + descriptif")
log_msg("=== ÉTAPES 8B.4-8B.5 TERMINÉES ===")
