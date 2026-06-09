# =============================================================================
# describe_concentration_monitoring.R
# Conc-1 — TDM monitoring coverage, selection bias, SMD
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(duckdb); library(DBI)
library(arrow)
library(ggplot2); library(patchwork)
library(lubridate)

log_msg("=== ÉTAPE 8B.1 : FAISABILITÉ DOSAGES VPA ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR, DIAGNOSTICS_DIR)

# Charger la cohorte principale
cohort <- read_parquet_dt(file.path(DATA_FINAL, "features_cohort.parquet"))
log_msg(sprintf("Cohorte dose-based: %d patients", nrow(cohort)))

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
dbExecute(con, "SET memory_limit='8GB'"); dbExecute(con, "SET threads=4")

lab_path  <- file.path(MIMIC_HOSP, "labevents.csv.gz")
VPA_ITEM  <- 51008L  # Valproic Acid — Blood — Chemistry

vpa_ids_str <- paste(cohort$subject_id, collapse=",")

# =============================================================================
# 1. Tous les dosages VPA pour la cohorte principale
# =============================================================================
log_msg("1. Extraction des dosages VPA de la cohorte...")

vpa_levels_q <- sprintf(
  "SELECT subject_id, itemid, charttime, valuenum, valueuom
   FROM read_csv_auto('%s', compression='gzip')
   WHERE subject_id IN (%s)
     AND itemid = %d
     AND valuenum IS NOT NULL AND valuenum > 0",
  lab_path, vpa_ids_str, VPA_ITEM
)
vpa_levels_raw <- as.data.table(dbGetQuery(con, vpa_levels_q))
vpa_levels_raw[, charttime := as.POSIXct(charttime, tz="UTC")]
log_msg(sprintf("Dosages VPA bruts: %d mesures pour %d patients",
                nrow(vpa_levels_raw), uniqueN(vpa_levels_raw$subject_id)))

# Vérification des unités
log_msg("Unités détectées:")
print(vpa_levels_raw[, .N, by=valueuom][order(-N)])

# Harmonisation : tout en µg/mL (µg/mL = mg/L numériquement identique)
vpa_levels_raw[, value_ugml := valuenum]  # mcg/mL déjà selon MIMIC convention

# =============================================================================
# 2. Fusionner avec time_zero de la cohorte
# =============================================================================
log_msg("2. Calcul du délai dosage → time zero...")

vpa_levels <- merge(vpa_levels_raw,
                     cohort[, .(subject_id, time_zero, high_dose, in_icu,
                                outcome_thrombo_100, plt_baseline)],
                     by="subject_id")
vpa_levels[, days_from_tz := as.numeric(difftime(charttime, time_zero, units="days"))]

# Statistiques de distribution temporelle
timing_summary <- vpa_levels[, .(
  n_measures     = .N,
  n_patients     = uniqueN(subject_id),
  days_min       = min(days_from_tz),
  days_p25       = quantile(days_from_tz, 0.25),
  days_median    = median(days_from_tz),
  days_p75       = quantile(days_from_tz, 0.75),
  days_max       = max(days_from_tz),
  pct_pre_tz     = mean(days_from_tz < 0)*100,
  pct_wk1        = mean(days_from_tz >= 0 & days_from_tz <= 7)*100,
  pct_wk2        = mean(days_from_tz > 7 & days_from_tz <= 14)*100,
  pct_wk4        = mean(days_from_tz > 14 & days_from_tz <= 30)*100
)]
log_msg("Distribution temporelle des dosages:")
print(timing_summary)

# Distribution fine
timing_by_window <- vpa_levels[, .(
  n_patients = uniqueN(subject_id),
  n_measures = .N
), by=.(window = cut(days_from_tz,
  breaks=c(-Inf,-7,-3,-1,0,1,2,3,5,7,10,14,21,30,Inf),
  labels=c("<-7j","-7:-3j","-3:-1j","-1:0j","J0-J1","J1-J2","J2-J3","J3-J5",
            "J5-J7","J7-J10","J10-J14","J14-J21","J21-J30",">J30")
))]
print(timing_by_window[order(window)])

# =============================================================================
# 3. Patients monitored vs non-monitored — bias assessment
# =============================================================================
log_msg("3. Comparison monitored vs non-monitored patients...")

# Patients monitorés = au moins 1 dosage dans les 7 premiers jours (J0-J7)
pts_monitored_early <- unique(vpa_levels[days_from_tz >= 0 & days_from_tz <= 7, subject_id])
cohort[, monitored_early := as.integer(subject_id %in% pts_monitored_early)]

# Aussi monitored = any level dans J0-J30
pts_monitored_30 <- unique(vpa_levels[days_from_tz >= 0 & days_from_tz <= 30, subject_id])
cohort[, monitored_30 := as.integer(subject_id %in% pts_monitored_30)]

# Aussi : dosage avant time zero (monitoring antérieur)
pts_monitored_pre <- unique(vpa_levels[days_from_tz < 0, subject_id])
cohort[, monitored_pre := as.integer(subject_id %in% pts_monitored_pre)]

log_msg(sprintf("Patients avec dosage J0-J7: %d (%.1f%%)",
                sum(cohort$monitored_early), mean(cohort$monitored_early)*100))
log_msg(sprintf("Patients avec dosage J0-J30: %d (%.1f%%)",
                sum(cohort$monitored_30), mean(cohort$monitored_30)*100))
log_msg(sprintf("Patients avec dosage pré-VPA: %d (%.1f%%)",
                sum(cohort$monitored_pre), mean(cohort$monitored_pre)*100))

# Caractéristiques monitored vs non-monitored
monitoring_compare <- cohort[, .(
  n              = .N,
  pct_highdose   = mean(high_dose)*100,
  pct_icu        = mean(in_icu, na.rm=T)*100,
  pct_sepsis     = mean(sepsis, na.rm=T)*100,
  pct_epilepsy   = mean(epilepsy_dx, na.rm=T)*100,
  pct_bipolar    = mean(bipolar_dx, na.rm=T)*100,
  pct_vpa_iv     = mean(vpa_iv, na.rm=T)*100,
  age_median     = median(age, na.rm=T),
  plt_baseline_med = median(plt_baseline, na.rm=T),
  alb_baseline_med = median(alb_baseline, na.rm=T),
  crea_baseline_med = median(crea_baseline, na.rm=T),
  daily_dose_med = median(daily_dose_main, na.rm=T),
  pct_thrombo    = mean(outcome_thrombo_100, na.rm=T)*100
), by=monitored_early]
print(monitoring_compare)
safe_write_csv(monitoring_compare, file.path(TABLES_DIR, "vpa_tdm_monitoring_vs_not.csv"), overwrite=TRUE)

# SMD monitored vs non-monitored
smd_vars <- c("high_dose","in_icu","sepsis","epilepsy_dx","vpa_iv","age","plt_baseline",
               "alb_baseline","crea_baseline","daily_dose_main","outcome_thrombo_100")
smd_monitoring <- sapply(intersect(smd_vars, names(cohort)), function(v) {
  tryCatch(smd(as.numeric(cohort[[v]]), cohort$monitored_early), error=function(e) NA_real_)
})
smd_mon_dt <- data.table(variable=names(smd_monitoring), smd=round(smd_monitoring, 3))
smd_mon_dt[, interpretation := fcase(
  abs(smd) < 0.1, "Balance OK",
  abs(smd) < 0.2, "Légère imbalance",
  default = "IMBALANCE SIGNIFICATIVE"
)]
log_msg("SMD monitored vs non-monitored:")
print(smd_mon_dt)
safe_write_csv(smd_mon_dt, file.path(TABLES_DIR, "vpa_tdm_monitoring_smd.csv"), overwrite=TRUE)

# =============================================================================
# 4. Distribution des concentrations
# =============================================================================
log_msg("4. Distribution des concentrations VPA...")

# Dans la fenêtre J0-J7 (première semaine post-initiation)
vpa_early <- vpa_levels[days_from_tz >= 0 & days_from_tz <= 7]
conc_summary <- vpa_early[, .(
  n_patients = uniqueN(subject_id),
  n_measures = .N,
  conc_min   = min(value_ugml),
  conc_p10   = quantile(value_ugml, 0.10),
  conc_p25   = quantile(value_ugml, 0.25),
  conc_median = median(value_ugml),
  conc_p75   = quantile(value_ugml, 0.75),
  conc_p90   = quantile(value_ugml, 0.90),
  conc_max   = max(value_ugml),
  pct_below50 = mean(value_ugml < 50)*100,
  pct_50_75  = mean(value_ugml >= 50 & value_ugml < 75)*100,
  pct_75_100 = mean(value_ugml >= 75 & value_ugml < 100)*100,
  pct_above100 = mean(value_ugml >= 100)*100
)]
log_msg("Distribution concentrations J0-J7:")
print(t(conc_summary))

safe_write_csv(conc_summary, file.path(TABLES_DIR, "vpa_tdm_feasibility.csv"), overwrite=TRUE)

# =============================================================================
# 5. Disponibilité conjointe dosage + outcome
# =============================================================================
log_msg("5. Disponibilité conjointe dosage + issues...")

# Pour les patients avec dosage J0-J7, vérifier co-disponibilité des labs d'outcome
# Charger les autres labs pour les patients monitored
all_item_ids_str <- paste(c(51265L, 50861L, 50878L, 50885L, 50862L, 50912L, 50866L), collapse=",")
monitored_ids_str <- paste(pts_monitored_early, collapse=",")

if (length(pts_monitored_early) > 0) {
  joint_q <- sprintf(
    "SELECT subject_id, itemid, COUNT(*) as n_measures,
            MIN(CASE WHEN itemid=51265 THEN valuenum END) as plt_min,
            MAX(CASE WHEN itemid=50861 THEN valuenum END) as alt_max,
            MAX(CASE WHEN itemid=50878 THEN valuenum END) as ast_max,
            MAX(CASE WHEN itemid=50866 THEN valuenum END) as nh3_max
     FROM read_csv_auto('%s', compression='gzip')
     WHERE subject_id IN (%s)
       AND itemid IN (%s)
       AND valuenum IS NOT NULL
     GROUP BY subject_id, itemid",
    lab_path, monitored_ids_str, all_item_ids_str
  )
  joint_avail <- as.data.table(dbGetQuery(con, joint_q))

  joint_by_item <- joint_avail[, .(n_patients=uniqueN(subject_id)), by=itemid]
  joint_by_item[, label := fcase(
    itemid==51265, "Plaquettes",
    itemid==50861, "ALT",
    itemid==50878, "AST",
    itemid==50885, "Bilirubin",
    itemid==50862, "Albumine",
    itemid==50912, "Créatinine",
    itemid==50866, "NH3"
  )]
  joint_by_item[, pct_monitored := n_patients / length(pts_monitored_early) * 100]
  log_msg("Co-disponibilité labs dans sous-cohorte monitorée:")
  print(joint_by_item)
  safe_write_csv(joint_by_item, file.path(TABLES_DIR, "vpa_tdm_joint_availability.csv"), overwrite=TRUE)
}

dbDisconnect(con, shutdown=TRUE)

# =============================================================================
# 6. Figures
# =============================================================================
log_msg("6. Figures...")

# Fig 1: Timing des dosages
p_timing <- ggplot(vpa_levels[days_from_tz >= -14 & days_from_tz <= 30],
                   aes(x=days_from_tz)) +
  geom_histogram(binwidth=0.5, fill="steelblue", alpha=0.8, color="white") +
  geom_vline(xintercept=0, color="red", linetype="dashed", linewidth=1) +
  annotate("rect", xmin=0, xmax=7, ymin=0, ymax=Inf,
           alpha=0.1, fill="green") +
  labs(title="Timing des dosages sanguins de valproate par rapport au time zero",
       subtitle="Ligne rouge = time zero (1ère prescription VPA) | Zone verte = fenêtre J0-J7",
       x="Jours depuis le time zero", y="N dosages") +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "vpa_tdm_timing.png"), p_timing, width=10, height=6, dpi=150)

# Fig 2: Distribution des concentrations
p_conc <- ggplot(vpa_early, aes(x=value_ugml)) +
  geom_histogram(bins=40, fill="steelblue", alpha=0.8, color="white") +
  geom_vline(xintercept=50, color="orange", linetype="dashed", linewidth=0.8) +
  geom_vline(xintercept=75, color="darkorange", linetype="dashed", linewidth=0.8) +
  geom_vline(xintercept=100, color="red", linetype="dashed", linewidth=0.8) +
  annotate("text", x=c(50,75,100), y=Inf, label=c("50","75","100 µg/mL"),
           angle=90, vjust=1.5, hjust=1.2, size=3.5, color=c("orange","darkorange","red")) +
  labs(title="Distribution des concentrations sanguines de valproate (J0-J7)",
       subtitle=sprintf("N=%d dosages chez %d patients | Lignes = seuils cliniques",
                        nrow(vpa_early), uniqueN(vpa_early$subject_id)),
       x="Concentration VPA (µg/mL)", y="N dosages") +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "concentration_distribution.png"), p_conc, width=10, height=6, dpi=150)

# Fig 3: SMD monitoring bias
p_smd_mon <- ggplot(smd_mon_dt[!is.na(smd)],
                     aes(x=smd, y=reorder(variable, abs(smd)),
                         color=interpretation)) +
  geom_point(size=3) +
  geom_vline(xintercept=c(-0.1, 0.1), color="gray50", linetype="dashed") +
  geom_vline(xintercept=0, color="black") +
  scale_color_manual(values=c("Balance OK"="steelblue",
                               "Légère imbalance"="orange",
                               "IMBALANCE SIGNIFICATIVE"="red")) +
  labs(title="Biais de monitoring : SMD monitored vs non-monitored (J0-J7)",
       subtitle="Patients avec dosage VPA vs sans dosage — sélection non aléatoire",
       x="Standardized Mean Difference", y="Variable", color="Interprétation") +
  theme_bw(base_size=12) + theme(legend.position="bottom")
ggsave(file.path(FIGURES_DIR, "vpa_tdm_monitoring_bias.png"),
       p_smd_mon, width=10, height=7, dpi=150)

# Fig 4: Concentration vs dose
conc_dose <- merge(
  vpa_early[, .(first_conc = value_ugml[which.min(days_from_tz)][1]), by=subject_id],
  cohort[, .(subject_id, daily_dose_main, high_dose, in_icu)],
  by="subject_id"
)
p_conc_dose <- ggplot(conc_dose[!is.na(first_conc) & first_conc < 200],
                       aes(x=daily_dose_main, y=first_conc, color=factor(high_dose))) +
  geom_point(alpha=0.4, size=1.5) +
  geom_smooth(method="loess", se=TRUE, color="black", linewidth=0.8) +
  geom_hline(yintercept=75, color="red", linetype="dashed") +
  geom_vline(xintercept=1000, color="blue", linetype="dashed") +
  scale_color_manual(values=c("steelblue","tomato"), labels=c("Faible dose","Haute dose")) +
  labs(title="Relation dose prescrite vs concentration mesurée (premier dosage J0-J7)",
       subtitle="Ligne rouge = seuil concentration 75 µg/mL | Ligne bleue = seuil dose 1000 mg/j",
       x="Dose journalière prescrite (mg/j)", y="Première concentration (µg/mL)",
       color="Groupe dose") +
  theme_bw(base_size=12)
ggsave(file.path(FIGURES_DIR, "concentration_vs_dose.png"), p_conc_dose, width=10, height=7, dpi=150)

timer_end(t0, "8B.1 faisabilité")

# Sauvegarder la cohorte enrichie avec flag monitoring
arrow::write_parquet(as.data.frame(cohort),
                      file.path(DATA_INTERMEDIATE, "cohort_with_monitoring_flag.parquet"))
arrow::write_parquet(as.data.frame(vpa_levels),
                      file.path(DATA_INTERMEDIATE, "vpa_levels_all.parquet"))
arrow::write_parquet(as.data.frame(vpa_early),
                      file.path(DATA_INTERMEDIATE, "vpa_levels_wk1.parquet"))

# =============================================================================
# Rapport 8B.1
# =============================================================================
n_mon  <- sum(cohort$monitored_early)
pct_mon <- mean(cohort$monitored_early)*100
n_tot  <- nrow(cohort)

report_8b1 <- sprintf(
'# 8B.1 Faisabilité Dosages VPA — Analyse Concentration-Based
**Date :** %s

## Identification du lab item

- **item_id MIMIC :** 51008 (`Valproic Acid`, Blood, Chemistry)
- **Unité :** µg/mL (identique à mg/L numériquement)
- **Plage de valeurs :** 3 – 510 µg/mL dans la cohorte VPA

## Disponibilité globale

- Patients avec au moins 1 dosage J0-J7 : **%d / %d (%.1f%%)**
- Patients avec au moins 1 dosage J0-J30 : **%d / %d (%.1f%%)**
- Mesures totales J0-J7 dans la cohorte : **%d**
- Médiane des concentrations J0-J7 : **%.0f µg/mL** [P25=%.0f, P75=%.0f]

## Timing des dosages

| Fenêtre | %% des mesures |
|---------|----------------|
| Avant time zero (pré-VPA) | %.1f%% |
| J0-J7 (1ère semaine) | %.1f%% |
| J7-J14 (2ème semaine) | %.1f%% |
| J14-J30 (3ème et 4ème semaines) | %.1f%% |

## Distribution des concentrations J0-J7

| Percentile | Valeur (µg/mL) |
|------------|----------------|
| P10 | %.0f |
| P25 | %.0f |
| P50 (médiane) | %.0f |
| P75 | %.0f |
| P90 | %.0f |
| %% < 50 µg/mL (sous-thérapeutique) | %.1f%% |
| %% 50-75 µg/mL (thérapeutique bas) | %.1f%% |
| %% 75-100 µg/mL (thérapeutique haut) | %.1f%% |
| %% > 100 µg/mL (suprathérapeutique) | %.1f%% |

## Biais de sélection des patients monitorés

Les patients avec dosage J0-J7 sont **différents** des non-monitorés :
%s

## Risques méthodologiques majeurs

1. **Biais de surveillance différentielle** : le dosage est prescrit en cas de suspicion de toxicité, d\'inefficacité, ou de changement de dose — sélection hautement informative
2. **Confounding by indication du TDM** : les patients les plus sévères, avec polymédication ou insuffisance organique, sont plus souvent dosés
3. **Timing inconnu trough/peak** : MIMIC ne distingue pas les dosages trough (creux) des dosages peak (pic). Le timing par rapport à la dernière administration est inférable mais imprécis
4. **Sélection sur survivants** : un patient doit rester hospitalisé et à risque jusqu\'au dosage — immortal time si non géré

## Recommandation de faisabilité

Avec **%.1f%%** des patients ayant un dosage exploitable en J0-J7, une analyse concentration-based est **faisable** dans une sous-cohorte.

Conception retenue :
- **Sous-cohorte monitorée** : patients avec dosage J0-J7
- **Exposition principale** : première concentration mesurée dans J0-J7
- **Seuil principal** : ≥ 75 µg/mL (limite supérieure de la zone thérapeutique conventionnelle basse)
- **Analyse de sensibilité** : ≥ 50 µg/mL, ≥ 100 µg/mL, quartiles, continue par splines

## Limites inhérentes

Le fait que seulement 60%% des patients soient monitorés — et que ces patients diffèrent systématiquement — rend toute inférence causale dans cette sous-cohorte sujette à caution.
Les résultats doivent être présentés comme **exploratoires** et **non directement comparables** à l\'analyse dose-based principale.
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  n_mon, n_tot, pct_mon,
  sum(cohort$monitored_30), n_tot, mean(cohort$monitored_30)*100,
  nrow(vpa_early),
  conc_summary$conc_median, conc_summary$conc_p25, conc_summary$conc_p75,
  timing_summary$pct_pre_tz, timing_summary$pct_wk1,
  timing_summary$pct_wk2, timing_summary$pct_wk4,
  conc_summary$conc_p10, conc_summary$conc_p25, conc_summary$conc_median,
  conc_summary$conc_p75, conc_summary$conc_p90,
  conc_summary$pct_below50, conc_summary$pct_50_75,
  conc_summary$conc_p90, conc_summary$pct_above100,
  {
    smd_sig <- smd_mon_dt[abs(smd) > 0.1 & !is.na(smd)]
    if (nrow(smd_sig) > 0)
      paste(sprintf("- `%s` : SMD = %.2f (%s)", smd_sig$variable, smd_sig$smd, smd_sig$interpretation),
            collapse="\n")
    else "Aucune imbalance majeure détectée"
  },
  pct_mon
)

writeLines(report_8b1, file.path(REPORTS_DIR, "08B_concentration_feasibility.md"))
log_msg("Rapport 8B.1 écrit")
log_msg("=== ÉTAPE 8B.1 TERMINÉE ===")
