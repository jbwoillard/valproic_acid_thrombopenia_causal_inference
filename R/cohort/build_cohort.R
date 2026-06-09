# =============================================================================
# build_cohort.R
# Step 1 — Cohort construction from MIMIC-IV
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(duckdb)
library(DBI)
library(arrow)
library(lubridate)

log_msg("=== ÉTAPE 5 : CONSTRUCTION DE LA COHORTE ===")
t0 <- timer_start()

ensure_dirs(DATA_FINAL, TABLES_DIR, REPORTS_DIR, DATA_INTERMEDIATE)

# IDs MIMIC-IV validés empiriquement dans l'étape de faisabilité
# Utiliser des IDs spécifiques (Blood uniquement) — confirmés dans d_labitems
platelets_id  <- 51265L   # Platelet Count — Blood — Hematology
alt_id        <- 50861L   # Alanine Aminotransferase (ALT) — Blood — Chemistry
ast_id        <- 50878L   # Asparate Aminotransferase (AST) — Blood — Chemistry
bilirubin_id  <- 50885L   # Bilirubin, Total — Blood — Chemistry
albumin_id    <- 50862L   # Albumin — Blood — Chemistry
creatinine_id <- 50912L   # Creatinine — Blood — Chemistry
ammonia_id    <- 50866L   # Ammonia — Blood — Chemistry
vpa_level_id  <- 51008L   # Valproic Acid — Blood — Chemistry

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "SET memory_limit='10GB'")
dbExecute(con, "SET threads=4")

presc_path  <- file.path(MIMIC_HOSP, "prescriptions.csv.gz")
lab_path    <- file.path(MIMIC_HOSP, "labevents.csv.gz")
adm_path    <- file.path(MIMIC_HOSP, "admissions.csv.gz")
pat_path    <- file.path(MIMIC_HOSP, "patients.csv.gz")
icu_path    <- file.path(MIMIC_ICU,  "icustays.csv.gz")
diag_path   <- file.path(MIMIC_HOSP, "diagnoses_icd.csv.gz")
svc_path    <- file.path(MIMIC_HOSP, "services.csv.gz")
omr_path    <- file.path(MIMIC_HOSP, "omr.csv.gz")
chart_path  <- file.path(MIMIC_ICU,  "chartevents.csv.gz")
input_path  <- file.path(MIMIC_ICU,  "inputevents.csv.gz")

# =============================================================================
# ÉTAPE A : Toutes les prescriptions VPA
# =============================================================================
log_msg("A. Extraction toutes prescriptions VPA...")

vpa_all_q <- sprintf(
  "SELECT subject_id, hadm_id, pharmacy_id,
          starttime, stoptime,
          drug, drug_type,
          CAST(dose_val_rx AS DOUBLE) as dose_val_rx,
          dose_unit_rx,
          CAST(doses_per_24_hrs AS DOUBLE) as doses_per_24_hrs,
          form_val_disp, form_unit_disp, route
   FROM read_csv_auto('%s', compression='gzip')
   WHERE (LOWER(drug) LIKE '%%valproat%%'
       OR LOWER(drug) LIKE '%%depakot%%'
       OR LOWER(drug) LIKE '%%depakene%%'
       OR LOWER(drug) LIKE '%%divalproex%%')
     AND drug_type = 'MAIN'",
  presc_path
)
vpa_all <- as.data.table(dbGetQuery(con, vpa_all_q))
vpa_all[, starttime := as.POSIXct(starttime, tz="UTC")]
vpa_all[, stoptime  := as.POSIXct(stoptime,  tz="UTC")]

n_step <- data.table(step = character(), n_admissions = integer(), n_patients = integer())

n_step <- rbind(n_step, data.table(
  step = "0. Prescriptions VPA brutes (MAIN)",
  n_admissions = uniqueN(vpa_all$hadm_id),
  n_patients   = uniqueN(vpa_all$subject_id)
))
log_msg(sprintf("A. VPA prescriptions: %d admissions, %d patients", n_step[1]$n_admissions, n_step[1]$n_patients))

# =============================================================================
# ÉTAPE B : Données patients et admissions
# =============================================================================
log_msg("B. Chargement patients et admissions...")

patients_q <- sprintf(
  "SELECT subject_id, gender, anchor_age, anchor_year, anchor_year_group, dod
   FROM read_csv_auto('%s', compression='gzip')",
  pat_path
)
patients <- as.data.table(dbGetQuery(con, patients_q))

adm_q <- sprintf(
  "SELECT subject_id, hadm_id, admittime, dischtime, deathtime,
          admission_type, admission_location, discharge_location,
          insurance, language, marital_status, race,
          hospital_expire_flag
   FROM read_csv_auto('%s', compression='gzip')",
  adm_path
)
admissions <- as.data.table(dbGetQuery(con, adm_q))
admissions[, admittime := as.POSIXct(admittime, tz="UTC")]
admissions[, dischtime := as.POSIXct(dischtime, tz="UTC")]

# Fusionner pour avoir l'âge et les infos demo
vpa_all <- merge(vpa_all, patients[, .(subject_id, gender, anchor_age)], by="subject_id", all.x=TRUE)
vpa_all <- merge(vpa_all, admissions[, .(hadm_id, admittime, dischtime, admission_type,
                                          hospital_expire_flag, race)],
                 by="hadm_id", all.x=TRUE)

# =============================================================================
# ÉTAPE C : Identifier la PREMIÈRE prescription VPA par admission (time zero)
# =============================================================================
log_msg("C. Calcul du time zero...")

time_zero <- vpa_all[
  !is.na(starttime),
  .(time_zero = min(starttime, na.rm=TRUE)),
  by = .(subject_id, hadm_id)
]

# Fusionner avec données admission
time_zero <- merge(time_zero, admissions[, .(hadm_id, admittime, dischtime, admission_type,
                                               hospital_expire_flag, race)],
                   by="hadm_id", all.x=TRUE)
time_zero <- merge(time_zero, patients[, .(subject_id, gender, anchor_age)], by="subject_id", all.x=TRUE)

n_step <- rbind(n_step, data.table(
  step = "1. Admissions avec time zero valide",
  n_admissions = nrow(time_zero),
  n_patients   = uniqueN(time_zero$subject_id)
))

# =============================================================================
# ÉTAPE D : Exclusion — âge < 18 ans
# =============================================================================
log_msg("D. Exclusion âge < 18...")
time_zero <- time_zero[anchor_age >= 18]
n_step <- rbind(n_step, data.table(
  step = "2. Exclure âge < 18 ans",
  n_admissions = nrow(time_zero),
  n_patients   = uniqueN(time_zero$subject_id)
))

# =============================================================================
# ÉTAPE E : Exclusion — prevalent users (VPA dans admission antérieure)
# =============================================================================
log_msg("E. Exclusion prevalent users...")

# Trier les admissions de chaque patient par time_zero
time_zero_sorted <- time_zero[order(subject_id, time_zero)]
# Identifier le rang de chaque admission par patient
time_zero_sorted[, adm_rank := seq_len(.N), by = subject_id]
# Ne garder que la PREMIÈRE admission avec VPA par patient (incident users)
incident_admissions <- time_zero_sorted[adm_rank == 1]

n_step <- rbind(n_step, data.table(
  step = "3. Exclure prevalent users (garder 1ère admission VPA par patient)",
  n_admissions = nrow(incident_admissions),
  n_patients   = uniqueN(incident_admissions$subject_id)
))
log_msg(sprintf("E. Après exclusion prevalent users: %d admissions/patients uniques",
                nrow(incident_admissions)))

# =============================================================================
# ÉTAPE F : Exclusion — séjour trop court (< 48h après time zero)
# =============================================================================
log_msg("F. Exclusion séjour < 48h après time zero...")
incident_admissions[, followup_h := as.numeric(difftime(dischtime, time_zero, units="hours"))]
incident_admissions_48h <- incident_admissions[followup_h >= 48 | is.na(followup_h)]
n_step <- rbind(n_step, data.table(
  step = "4. Exclure séjour < 48h après time zero",
  n_admissions = nrow(incident_admissions_48h),
  n_patients   = uniqueN(incident_admissions_48h$subject_id)
))

# =============================================================================
# ÉTAPE G : Calcul de la dose journalière (72h après time zero)
# =============================================================================
log_msg("G. Calcul dose journalière (fenêtre 72h)...")

# Joindre les prescriptions à la cohorte
vpa_dose <- merge(vpa_all[, .(subject_id, hadm_id, starttime, stoptime,
                                dose_val_rx, dose_unit_rx, doses_per_24_hrs, route)],
                  incident_admissions_48h[, .(subject_id, hadm_id, time_zero)],
                  by=c("subject_id","hadm_id"))

# Filtrer prescriptions dans la fenêtre 0-72h après time zero
vpa_dose[, hours_from_tz := as.numeric(difftime(starttime, time_zero, units="hours"))]
vpa_dose_window <- vpa_dose[hours_from_tz >= 0 & hours_from_tz <= 72 & !is.na(dose_val_rx) & dose_val_rx > 0]

# Convertir unités (g → mg si nécessaire)
vpa_dose_window[dose_unit_rx == "g", dose_val_rx := dose_val_rx * 1000]
vpa_dose_window[dose_unit_rx == "g", dose_unit_rx := "mg"]

# Calculer la dose journalière par prescription
vpa_dose_window[, dose_per_day := dose_val_rx * ifelse(!is.na(doses_per_24_hrs) & doses_per_24_hrs > 0,
                                                        doses_per_24_hrs, 1)]

# Agréger par admission : dose journalière maximale, moyenne, et somme totale dans la fenêtre
dose_summary <- vpa_dose_window[, .(
  daily_dose_max     = max(dose_per_day, na.rm=TRUE),
  daily_dose_mean    = mean(dose_per_day, na.rm=TRUE),
  total_dose_72h     = sum(dose_val_rx, na.rm=TRUE),
  n_prescriptions_72h = .N,
  has_iv             = any(route %in% c("IV", "IV DRIP", "IVPB"), na.rm=TRUE)
), by = .(subject_id, hadm_id)]

dose_summary[, daily_dose_main := daily_dose_max]
dose_summary[, high_dose := as.integer(daily_dose_main > VPA_HIGH_DOSE_THRESHOLD_MG)]

# Joindre à la cohorte
cohort <- merge(incident_admissions_48h, dose_summary, by=c("subject_id","hadm_id"), all.x=TRUE)

# Exclure ceux sans dose calculable dans la fenêtre 72h
cohort_dose_ok <- cohort[!is.na(daily_dose_main) & daily_dose_main > 0]
n_step <- rbind(n_step, data.table(
  step = "5. Avec dose calculable dans les 72h après time zero",
  n_admissions = nrow(cohort_dose_ok),
  n_patients   = uniqueN(cohort_dose_ok$subject_id)
))
log_msg(sprintf("G. Après calcul dose: %d admissions, %d patients", nrow(cohort_dose_ok), uniqueN(cohort_dose_ok$subject_id)))

# =============================================================================
# ÉTAPE H : Extraction des labs baseline ([-7j, 0j] avant time zero)
# =============================================================================
log_msg("H. Extraction labs baseline...")

vpa_ids_str <- paste(unique(cohort_dose_ok$subject_id), collapse=",")
all_item_ids_str <- paste(c(platelets_id, alt_id, ast_id, bilirubin_id, albumin_id,
                              creatinine_id, ammonia_id, vpa_level_id), collapse=",")

# Récupérer tous les labs de la cohorte (par subject_id seulement, hadm_id peut être NULL)
labs_q <- sprintf(
  "SELECT subject_id, itemid, charttime,
          valuenum, valueuom
   FROM read_csv_auto('%s', compression='gzip')
   WHERE subject_id IN (%s)
     AND itemid IN (%s)
     AND valuenum IS NOT NULL
     AND valuenum > 0",
  lab_path, vpa_ids_str, all_item_ids_str
)
labs_all <- as.data.table(dbGetQuery(con, labs_q))
labs_all[, charttime := as.POSIXct(charttime, tz="UTC")]
log_msg(sprintf("H. Labs récupérés: %d mesures pour %d patients", nrow(labs_all), uniqueN(labs_all$subject_id)))

# Fusionner avec time_zero par subject_id (un time_zero par patient dans cette cohorte)
labs_all <- merge(labs_all, cohort_dose_ok[, .(subject_id, hadm_id, time_zero)],
                  by="subject_id")
labs_all[, days_from_tz := as.numeric(difftime(charttime, time_zero, units="days"))]

# Fonction générique : dernière valeur baseline dans [-window_days, 0)
get_baseline_value <- function(dt, item_ids, window_days = 7, col_name = "value") {
  baseline <- dt[
    itemid %in% item_ids &
    days_from_tz >= -window_days &
    days_from_tz < 0
  ]
  if (nrow(baseline) == 0) return(data.table(subject_id = integer()))
  baseline[, .SD[which.max(days_from_tz)], by=.(subject_id)][,
    .(subject_id, value = valuenum)
  ] |> setnames("value", col_name)
}

bl_plt  <- get_baseline_value(labs_all, platelets_id,  7, "plt_baseline")
bl_alt  <- get_baseline_value(labs_all, alt_id,        7, "alt_baseline")
bl_ast  <- get_baseline_value(labs_all, ast_id,        7, "ast_baseline")
bl_bili <- get_baseline_value(labs_all, bilirubin_id,  7, "bili_baseline")
bl_alb  <- get_baseline_value(labs_all, albumin_id,    7, "alb_baseline")
bl_crea <- get_baseline_value(labs_all, creatinine_id, 7, "crea_baseline")
bl_nh3  <- get_baseline_value(labs_all, ammonia_id,    7, "nh3_baseline")
bl_vpa  <- get_baseline_value(labs_all, vpa_level_id,  7, "vpa_level_baseline")

for (bl in list(bl_plt, bl_alt, bl_ast, bl_bili, bl_alb, bl_crea, bl_nh3, bl_vpa)) {
  if (nrow(bl) > 0) {
    cohort_dose_ok <- merge(cohort_dose_ok, bl, by="subject_id", all.x=TRUE)
  }
}

log_msg(sprintf("H. Patients avec plaquettes baseline: %.0f (%.0f%%)",
                sum(!is.na(cohort_dose_ok$plt_baseline)),
                round(mean(!is.na(cohort_dose_ok$plt_baseline))*100)))

# =============================================================================
# ÉTAPE I : Exclusion — pas de plaquettes baseline (outcome principal non évaluable)
# =============================================================================
log_msg("I. Exclusion sans plaquettes baseline (fenêtre 7j)...")
cohort_plt_bl <- cohort_dose_ok[!is.na(plt_baseline)]
n_step <- rbind(n_step, data.table(
  step = "6. Avec plaquettes baseline disponibles (J-7 à J0)",
  n_admissions = nrow(cohort_plt_bl),
  n_patients   = uniqueN(cohort_plt_bl$subject_id)
))

# =============================================================================
# ÉTAPE J : Exclusion — thrombopénie préexistante (pour outcome principal)
# =============================================================================
log_msg("J. Exclusion thrombopénie préexistante (plt_baseline < 100)...")
cohort_no_prior_thrombo <- cohort_plt_bl[plt_baseline >= 100]
n_step <- rbind(n_step, data.table(
  step = "7. Exclure thrombopénie préexistante (plt baseline < 100 G/L)",
  n_admissions = nrow(cohort_no_prior_thrombo),
  n_patients   = uniqueN(cohort_no_prior_thrombo$subject_id)
))
log_msg(sprintf("J. Cohorte principale (thrombo): %d admissions, %d patients",
                nrow(cohort_no_prior_thrombo), uniqueN(cohort_no_prior_thrombo$subject_id)))

# =============================================================================
# ÉTAPE K : Définition des outcomes
# =============================================================================
log_msg("K. Définition des outcomes (30j post time zero)...")

# Outcome principal : thrombopénie
outcome_plt <- labs_all[
  itemid == platelets_id &
  days_from_tz >= 0 &
  days_from_tz <= MAX_FOLLOWUP_DAYS
]
# Plus mauvaise valeur par patient dans le suivi
min_plt_fu <- outcome_plt[, .(
  plt_min_fu  = min(valuenum, na.rm=TRUE),
  plt_fu_time = min(days_from_tz[valuenum == min(valuenum, na.rm=TRUE)][1], na.rm=TRUE),
  plt_n_fu    = .N
), by=.(subject_id)]

cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, min_plt_fu, by="subject_id", all.x=TRUE)

# Indicateurs d'outcome
cohort_no_prior_thrombo[, outcome_thrombo_100 := as.integer(!is.na(plt_min_fu) & plt_min_fu < 100)]
cohort_no_prior_thrombo[, outcome_thrombo_50  := as.integer(!is.na(plt_min_fu) & plt_min_fu < 50)]
cohort_no_prior_thrombo[, outcome_thrombo_rel30 :=
  as.integer(!is.na(plt_min_fu) & !is.na(plt_baseline) &
             (plt_baseline - plt_min_fu) / plt_baseline > 0.30)]
cohort_no_prior_thrombo[, outcome_thrombo_rel50 :=
  as.integer(!is.na(plt_min_fu) & !is.na(plt_baseline) &
             (plt_baseline - plt_min_fu) / plt_baseline > 0.50)]

# Outcome ALT/AST (hépatotoxicité)
outcome_alt <- labs_all[itemid == alt_id & days_from_tz >= 0 & days_from_tz <= MAX_FOLLOWUP_DAYS,
  .(alt_max_fu = max(valuenum, na.rm=TRUE), alt_n_fu = .N), by=.(subject_id)]
outcome_ast <- labs_all[itemid == ast_id & days_from_tz >= 0 & days_from_tz <= MAX_FOLLOWUP_DAYS,
  .(ast_max_fu = max(valuenum, na.rm=TRUE)), by=.(subject_id)]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, outcome_alt, by="subject_id", all.x=TRUE)
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, outcome_ast, by="subject_id", all.x=TRUE)
cohort_no_prior_thrombo[, outcome_hepato_3xULN :=
  as.integer((!is.na(alt_max_fu) & alt_max_fu > 120) | (!is.na(ast_max_fu) & ast_max_fu > 120))]

# Outcome NH3 (hyperammoniémie)
outcome_nh3 <- labs_all[itemid == ammonia_id & days_from_tz >= 0 & days_from_tz <= MAX_FOLLOWUP_DAYS,
  .(nh3_max_fu = max(valuenum, na.rm=TRUE)), by=.(subject_id)]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, outcome_nh3, by="subject_id", all.x=TRUE)
cohort_no_prior_thrombo[, outcome_hypernh3 :=
  as.integer(!is.na(nh3_max_fu) & nh3_max_fu > 55)]

# Outcome composite exploratoire
cohort_no_prior_thrombo[, outcome_composite :=
  as.integer(outcome_thrombo_100 == 1 | outcome_hepato_3xULN == 1)]

log_msg(sprintf("K. Outcomes: thrombo100=%.0f%%, hepato=%.0f%%, NH3=%.1f%%",
                mean(cohort_no_prior_thrombo$outcome_thrombo_100, na.rm=T)*100,
                mean(cohort_no_prior_thrombo$outcome_hepato_3xULN, na.rm=T)*100,
                mean(cohort_no_prior_thrombo$outcome_hypernh3, na.rm=T)*100))

# =============================================================================
# ÉTAPE L : Variables supplémentaires (ICU, comorbidités)
# =============================================================================
log_msg("L. Variables ICU et comorbidités...")

# ICU
icu_q <- sprintf(
  "SELECT subject_id, hadm_id,
          MIN(intime) as icu_first_intime,
          MAX(los) as icu_max_los,
          COUNT(*) as n_icu_stays
   FROM read_csv_auto('%s', compression='gzip')
   WHERE subject_id IN (%s)
   GROUP BY subject_id, hadm_id",
  icu_path, vpa_ids_str
)
icu_stays <- as.data.table(dbGetQuery(con, icu_q))
icu_stays[, in_icu := 1L]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, icu_stays, by=c("subject_id","hadm_id"), all.x=TRUE)
cohort_no_prior_thrombo[is.na(in_icu), in_icu := 0L]

# Diagnoses — sepsis (ICD-10 A41*, ICD-9 038*)
diag_q <- sprintf(
  "SELECT DISTINCT subject_id, hadm_id, icd_code, icd_version
   FROM read_csv_auto('%s', compression='gzip')
   WHERE subject_id IN (%s)",
  diag_path, vpa_ids_str
)
diagnoses <- as.data.table(dbGetQuery(con, diag_q))

# Sepsis
sepsis_icd10 <- diagnoses[icd_version == 10 & grepl("^A41", icd_code), .(subject_id, hadm_id)][, sepsis := 1L]
sepsis_icd9  <- diagnoses[icd_version == 9  & grepl("^038", icd_code), .(subject_id, hadm_id)][, sepsis := 1L]
sepsis_all   <- unique(rbind(sepsis_icd10, sepsis_icd9))
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, sepsis_all, by=c("subject_id","hadm_id"), all.x=TRUE)
cohort_no_prior_thrombo[is.na(sepsis), sepsis := 0L]

# Epilepsie (ICD-10 G40*, ICD-9 345*)
epilepsy_icd <- diagnoses[
  (icd_version == 10 & grepl("^G40", icd_code)) |
  (icd_version == 9  & grepl("^345", icd_code)),
  .(subject_id, hadm_id)
][, epilepsy_dx := 1L]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, unique(epilepsy_icd), by=c("subject_id","hadm_id"), all.x=TRUE)
cohort_no_prior_thrombo[is.na(epilepsy_dx), epilepsy_dx := 0L]

# Trouble bipolaire (ICD-10 F31*, ICD-9 296.0-296.1)
bipolar_icd <- diagnoses[
  (icd_version == 10 & grepl("^F31", icd_code)) |
  (icd_version == 9  & grepl("^296[01]", icd_code)),
  .(subject_id, hadm_id)
][, bipolar_dx := 1L]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, unique(bipolar_icd), by=c("subject_id","hadm_id"), all.x=TRUE)
cohort_no_prior_thrombo[is.na(bipolar_dx), bipolar_dx := 0L]

# =============================================================================
# ÉTAPE M : Co-médications baseline (prescriptions avant time zero)
# =============================================================================
log_msg("M. Co-médications baseline...")

# Toutes prescriptions dans les 7j avant time zero
comed_q <- sprintf(
  "SELECT p.subject_id, p.hadm_id,
          LOWER(p.drug) as drug_lower
   FROM read_csv_auto('%s', compression='gzip') p
   INNER JOIN (SELECT hadm_id, time_zero FROM
     (VALUES %s) AS tz(hadm_id, time_zero)) t ON p.hadm_id = t.hadm_id
   WHERE p.starttime <= t.time_zero
     AND p.starttime >= CAST(t.time_zero AS TIMESTAMP) - INTERVAL '7 days'
     AND p.drug_type = 'MAIN'",
  presc_path,
  paste(sprintf("(%d, '%s')",
                cohort_no_prior_thrombo$hadm_id,
                format(cohort_no_prior_thrombo$time_zero, "%Y-%m-%d %H:%M:%S")),
        collapse=", ")
)

tryCatch({
  comed_raw <- as.data.table(dbGetQuery(con, comed_q))

  # Médicaments thrombopéniants
  thrombopenia_drugs <- c("heparin", "enoxaparin", "fondaparinux", "linezolid",
                           "vancomycin", "trimethoprim", "platelet")
  # Hépatotoxiques
  hepatotox_drugs <- c("acetaminophen", "amiodarone", "methotrexate", "isoniazid",
                        "fluconazole", "rifampin", "nitrofurantoin")
  # Interactions VPA
  carbapenem_drugs <- c("meropenem", "imipenem", "ertapenem", "doripenem")
  topiramate_drugs <- c("topiramate", "topamax")

  drug_flags <- comed_raw[, .(
    comed_topiramate   = as.integer(any(grepl(paste(topiramate_drugs, collapse="|"), drug_lower))),
    comed_carbapenem   = as.integer(any(grepl(paste(carbapenem_drugs, collapse="|"), drug_lower))),
    comed_heparin      = as.integer(any(grepl("heparin|enoxaparin|fondaparinux", drug_lower))),
    comed_linezolid    = as.integer(any(grepl("linezolid", drug_lower))),
    comed_steroid      = as.integer(any(grepl("methyl|prednis|dexameth|hydrocort", drug_lower))),
    comed_antipsy      = as.integer(any(grepl("quetiapine|olanzapine|risperidone|haloperidol|clozapine", drug_lower))),
    comed_other_aed    = as.integer(any(grepl("phenytoin|levetiracetam|lamotrigine|carbamazepine|phenobarbital", drug_lower)))
  ), by=.(subject_id, hadm_id)]

  cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, drug_flags, by=c("subject_id","hadm_id"), all.x=TRUE)
  for (col in c("comed_topiramate","comed_carbapenem","comed_heparin","comed_linezolid",
                "comed_steroid","comed_antipsy","comed_other_aed")) {
    cohort_no_prior_thrombo[is.na(get(col)), (col) := 0L]
  }
  log_msg("M. Co-médications extraites avec succès")
}, error = function(e) {
  log_msg(sprintf("M. Erreur co-médications (stratégie alternative): %s", e$message), level="WARN")
  # Si la requête VALUES est trop grande, utiliser une approche plus simple
  all_presc_q <- sprintf(
    "SELECT subject_id, hadm_id, LOWER(drug) as drug_lower
     FROM read_csv_auto('%s', compression='gzip')
     WHERE subject_id IN (%s) AND drug_type = 'MAIN'",
    presc_path, vpa_ids_str
  )
  comed_all <- as.data.table(dbGetQuery(con, all_presc_q))
  comed_all <- merge(comed_all,
                     cohort_no_prior_thrombo[, .(subject_id, hadm_id, time_zero)],
                     by=c("subject_id","hadm_id"))
  drug_flags <- comed_all[, .(
    comed_topiramate = as.integer(any(grepl("topiramate|topamax", drug_lower))),
    comed_carbapenem = as.integer(any(grepl("meropenem|imipenem|ertapenem|doripenem", drug_lower))),
    comed_heparin    = as.integer(any(grepl("heparin|enoxaparin|fondaparinux", drug_lower))),
    comed_linezolid  = as.integer(any(grepl("linezolid", drug_lower))),
    comed_steroid    = as.integer(any(grepl("methyl|prednis|dexameth|hydrocort", drug_lower))),
    comed_antipsy    = as.integer(any(grepl("quetiapine|olanzapine|risperidone|haloperidol", drug_lower))),
    comed_other_aed  = as.integer(any(grepl("phenytoin|levetiracetam|lamotrigine|carbamazepine|phenobarbital", drug_lower)))
  ), by=.(subject_id, hadm_id)]
  cohort_no_prior_thrombo <<- merge(cohort_no_prior_thrombo, drug_flags, by=c("subject_id","hadm_id"), all.x=TRUE)
  for (col in c("comed_topiramate","comed_carbapenem","comed_heparin","comed_linezolid",
                "comed_steroid","comed_antipsy","comed_other_aed")) {
    cohort_no_prior_thrombo[is.na(get(col)), (col) := 0L]
  }
  log_msg("M. Co-médications extraites (méthode alternative)")
})

# =============================================================================
# ÉTAPE N : Poids baseline (OMR)
# =============================================================================
log_msg("N. Poids baseline depuis OMR...")
weight_q <- sprintf(
  "SELECT subject_id, chartdate, result_value
   FROM read_csv_auto('%s', compression='gzip')
   WHERE subject_id IN (%s)
     AND result_name = 'Weight (Lbs)'
     AND result_value IS NOT NULL",
  omr_path, vpa_ids_str
)
weight_omr <- as.data.table(dbGetQuery(con, weight_q))
weight_omr[, chartdate := as.Date(chartdate)]
weight_omr[, weight_kg := as.numeric(result_value) * 0.453592]  # lbs → kg
# Fusionner avec cohorte (prendre le poids le plus proche du time zero)
weight_omr <- merge(weight_omr, cohort_no_prior_thrombo[, .(subject_id, time_zero)],
                    by="subject_id", all.x=TRUE)
weight_omr[, days_from_tz := as.numeric(as.Date(time_zero) - chartdate)]
weight_baseline_omr <- weight_omr[
  days_from_tz >= 0 & days_from_tz <= 180  # poids dans les 6 mois avant
][, .SD[which.min(days_from_tz)], by=subject_id][,
  .(subject_id, weight_kg_baseline = weight_kg)
]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, weight_baseline_omr, by="subject_id", all.x=TRUE)

# =============================================================================
# ÉTAPE O : Variables finales et nettoyage
# =============================================================================
log_msg("O. Finalisation de la cohorte...")

# Variables dérivées
cohort_no_prior_thrombo[, age := anchor_age]
cohort_no_prior_thrombo[, female := as.integer(gender == "F")]
cohort_no_prior_thrombo[, emergency := as.integer(grepl("EMERGENCY|URGENT", toupper(admission_type)))]
cohort_no_prior_thrombo[, followup_days := pmin(as.numeric(difftime(dischtime, time_zero, units="days")),
                                                 MAX_FOLLOWUP_DAYS, na.rm=TRUE)]
cohort_no_prior_thrombo[, in_hospital_death := hospital_expire_flag]

# Durée de préhospitalisation
cohort_no_prior_thrombo[, pre_hosp_days := as.numeric(difftime(time_zero, admittime, units="days"))]
cohort_no_prior_thrombo[pre_hosp_days < 0, pre_hosp_days := 0]

# Exposition secondaire : niveau sérique VPA dans la semaine suivant time zero
vpa_level_fu <- labs_all[
  itemid == vpa_level_id &
  days_from_tz >= 0 & days_from_tz <= 7,
  .(vpa_level_mean_wk1 = mean(valuenum, na.rm=TRUE),
    vpa_level_max_wk1  = max(valuenum, na.rm=TRUE),
    n_vpa_levels_wk1   = .N),
  by=.(subject_id)
]
cohort_no_prior_thrombo <- merge(cohort_no_prior_thrombo, vpa_level_fu, by="subject_id", all.x=TRUE)
cohort_no_prior_thrombo[, has_vpa_level := as.integer(!is.na(vpa_level_mean_wk1))]

# Tertiles de dose pour analyses HTE
cohort_no_prior_thrombo[, dose_tertile := cut(daily_dose_main,
  breaks = quantile(daily_dose_main, c(0, 1/3, 2/3, 1), na.rm=TRUE),
  labels = c("low", "medium", "high"),
  include.lowest = TRUE)]

dbDisconnect(con, shutdown=TRUE)

# =============================================================================
# ÉTAPE P : Flowchart final
# =============================================================================
log_msg("P. Flowchart final...")
n_step <- rbind(n_step, data.table(
  step = "8. Cohorte analytique finale (outcome principal : thrombopénie)",
  n_admissions = nrow(cohort_no_prior_thrombo),
  n_patients   = uniqueN(cohort_no_prior_thrombo$subject_id)
))

# Décomposition haute vs basse dose
n_high  <- sum(cohort_no_prior_thrombo$high_dose == 1, na.rm=TRUE)
n_low   <- sum(cohort_no_prior_thrombo$high_dose == 0, na.rm=TRUE)
log_msg(sprintf("P. Haute dose (>1000 mg/j): %d (%.1f%%) | Faible dose: %d (%.1f%%)",
                n_high, n_high/nrow(cohort_no_prior_thrombo)*100,
                n_low,  n_low/nrow(cohort_no_prior_thrombo)*100))
log_msg(sprintf("P. Outcome thrombopénie: %d (%.1f%%)",
                sum(cohort_no_prior_thrombo$outcome_thrombo_100, na.rm=TRUE),
                mean(cohort_no_prior_thrombo$outcome_thrombo_100, na.rm=TRUE)*100))

print(n_step)
safe_write_csv(n_step, file.path(TABLES_DIR, "cohort_flow.csv"), overwrite=TRUE)

# =============================================================================
# SAUVEGARDE
# =============================================================================
log_msg("Q. Sauvegarde cohorte finale...")
safe_write_parquet(cohort_no_prior_thrombo,
                   file.path(DATA_FINAL, "analytic_cohort.parquet"), overwrite=TRUE)
safe_write_csv(cohort_no_prior_thrombo,
               file.path(DATA_FINAL, "analytic_cohort.csv"), overwrite=TRUE)

# Rapport
timer_end(t0, "Étape 5 totale")

cohort_report <- sprintf(
'# Construction de la Cohorte Analytique — Étude VPA/MIMIC-IV
**Date :** %s

## Flowchart de Sélection

| Étape | N admissions | N patients |
|-------|-------------|------------|
%s

## Cohorte Finale

- **N patients :** %d
- **N admissions :** %d
- **Haute dose (> 1000 mg/j) :** %d (%.1f%%)
- **Faible dose (≤ 1000 mg/j) :** %d (%.1f%%)

## Outcomes

| Outcome | N | Prévalence |
|---------|---|-----------|
| Thrombopénie < 100 G/L | %d | %.1f%% |
| Thrombopénie < 50 G/L | %d | %.1f%% |
| Baisse plaquettes > 30%% | %d | %.1f%% |
| Hépatotoxicité (ALT/AST > 3×ULN) | %d | %.1f%% |
| Hyperammoniémie (NH3 > 55) | %d | %.1f%% |
| Composite (thrombo OR hepato) | %d | %.1f%% |

## Disponibilité des covariables

| Variable | N disponible | %% |
|----------|-------------|-----|
| Plaquettes baseline | %d | %.0f%% |
| ALT baseline | %d | %.0f%% |
| AST baseline | %d | %.0f%% |
| Albumine baseline | %d | %.0f%% |
| Créatinine baseline | %d | %.0f%% |
| Dosage VPA sérique (J0-J7) | %d | %.0f%% |
| Poids baseline | %d | %.0f%% |

## Notes méthodologiques

1. Design incident user : seule la première admission avec VPA par patient est retenue
2. Time zero = première prescription VPA (MAIN type)
3. La dose journalière est calculée sur une fenêtre de 72h après time zero
4. Les plaquettes baselines sont cherchées dans la fenêtre [-7j, 0j]
5. Le suivi est défini sur 30 jours ou jusqu'\''à la sortie si antérieure
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(apply(n_step, 1, function(r)
    sprintf("| %s | %s | %s |",
            r["step"],
            format(as.integer(r["n_admissions"]), big.mark=","),
            format(as.integer(r["n_patients"]), big.mark=","))),
    collapse="\n"),
  nrow(cohort_no_prior_thrombo),
  nrow(cohort_no_prior_thrombo),
  n_high, n_high/nrow(cohort_no_prior_thrombo)*100,
  n_low,  n_low/nrow(cohort_no_prior_thrombo)*100,
  sum(cohort_no_prior_thrombo$outcome_thrombo_100, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_thrombo_100, na.rm=T)*100,
  sum(cohort_no_prior_thrombo$outcome_thrombo_50, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_thrombo_50, na.rm=T)*100,
  sum(cohort_no_prior_thrombo$outcome_thrombo_rel30, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_thrombo_rel30, na.rm=T)*100,
  sum(cohort_no_prior_thrombo$outcome_hepato_3xULN, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_hepato_3xULN, na.rm=T)*100,
  sum(cohort_no_prior_thrombo$outcome_hypernh3, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_hypernh3, na.rm=T)*100,
  sum(cohort_no_prior_thrombo$outcome_composite, na.rm=T),
  mean(cohort_no_prior_thrombo$outcome_composite, na.rm=T)*100,
  sum(!is.na(cohort_no_prior_thrombo$plt_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$plt_baseline))*100,
  sum(!is.na(cohort_no_prior_thrombo$alt_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$alt_baseline))*100,
  sum(!is.na(cohort_no_prior_thrombo$ast_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$ast_baseline))*100,
  sum(!is.na(cohort_no_prior_thrombo$alb_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$alb_baseline))*100,
  sum(!is.na(cohort_no_prior_thrombo$crea_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$crea_baseline))*100,
  sum(cohort_no_prior_thrombo$has_vpa_level, na.rm=T),
  mean(cohort_no_prior_thrombo$has_vpa_level, na.rm=T)*100,
  sum(!is.na(cohort_no_prior_thrombo$weight_kg_baseline)),
  mean(!is.na(cohort_no_prior_thrombo$weight_kg_baseline))*100
)
writeLines(cohort_report, file.path(REPORTS_DIR, "05_cohort_construction.md"))
log_msg("Rapport cohorte écrit:", file.path(REPORTS_DIR, "05_cohort_construction.md"))
log_msg("=== ÉTAPE 5 TERMINÉE ===")
