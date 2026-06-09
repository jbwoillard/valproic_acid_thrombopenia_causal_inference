# =============================================================================
# engineer_features.R
# Step 2 — Feature engineering, MICE imputation, propensity variables
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(arrow)
library(mice)
library(naniar)

log_msg("=== ÉTAPE 6 : FEATURE ENGINEERING ===")
t0 <- timer_start()

# Charger la cohorte
cohort <- read_parquet_dt(file.path(DATA_FINAL, "analytic_cohort.parquet"))
log_msg(sprintf("Cohorte chargée: %d patients, %d variables", nrow(cohort), ncol(cohort)))

# =============================================================================
# 1. Nettoyage et encodage des variables
# =============================================================================
log_msg("1. Encodage des variables...")

# Âge
cohort[, age := as.numeric(anchor_age)]

# Sexe
cohort[, female := as.integer(gender == "F")]

# Race (simplifiée)
cohort[, race_simple := fcase(
  grepl("WHITE", toupper(race)), "white",
  grepl("BLACK|AFRICAN", toupper(race)), "black",
  grepl("ASIAN", toupper(race)), "asian",
  grepl("HISPANIC|LATIN", toupper(race)), "hispanic",
  default = "other_unknown"
)]

# Urgence
cohort[, emergency := as.integer(grepl("EMERGENCY|URGENT", toupper(admission_type)))]

# ICU (déjà binaire)
cohort[is.na(in_icu), in_icu := 0L]

# Durée pré-hospitalisation (jours avant VPA)
cohort[, pre_hosp_days_winsorized := pmin(pre_hosp_days, 30, na.rm=TRUE)]

# =============================================================================
# 2. Transformation des variables biologiques baselines
# =============================================================================
log_msg("2. Variables biologiques baseline...")

# Albumine (g/dL)
cohort[, alb_baseline_low := as.integer(!is.na(alb_baseline) & alb_baseline < 3.0)]
cohort[, alb_group := fcase(
  is.na(alb_baseline), "missing",
  alb_baseline < 2.5, "severe_low",
  alb_baseline < 3.0, "low",
  alb_baseline < 3.5, "borderline",
  default = "normal"
)]

# Créatinine (mg/dL) → insuffisance rénale
cohort[, crea_aki := as.integer(!is.na(crea_baseline) & crea_baseline >= 1.5)]
cohort[, crea_group := fcase(
  is.na(crea_baseline), "missing",
  crea_baseline < 1.0, "normal",
  crea_baseline < 1.5, "mild",
  crea_baseline < 2.0, "moderate",
  default = "severe"
)]

# Plaquettes baseline (G/L) — déjà filtrée >= 100
cohort[, plt_baseline_100_150 := as.integer(!is.na(plt_baseline) & plt_baseline >= 100 & plt_baseline < 150)]
cohort[, plt_baseline_group := fcase(
  is.na(plt_baseline), "missing",
  plt_baseline < 150, "low_normal",
  plt_baseline < 250, "normal",
  default = "high"
)]

# ALT/AST baseline — indicateur d'hépatopathie préexistante
cohort[, liver_baseline_elevated := as.integer(
  (!is.na(alt_baseline) & alt_baseline > 40) |
  (!is.na(ast_baseline) & ast_baseline > 40)
)]

# =============================================================================
# 3. Exposition : dose et variantes
# =============================================================================
log_msg("3. Variables d'exposition...")

# Dose principale (déjà calculée)
cohort[, log_dose := log(daily_dose_main + 1)]
cohort[, dose_per_kg := ifelse(!is.na(weight_kg_baseline) & weight_kg_baseline > 0,
                                daily_dose_main / weight_kg_baseline, NA_real_)]

# Variables catégorielles de dose
cohort[, dose_cat_750  := as.integer(daily_dose_main > 750)]
cohort[, dose_cat_1000 := as.integer(daily_dose_main > 1000)]  # principale
cohort[, dose_cat_1250 := as.integer(daily_dose_main > 1250)]
cohort[, dose_cat_1500 := as.integer(daily_dose_main > 1500)]

# Tertiles de dose (déjà dans la cohorte comme dose_tertile)
# Quartiles de dose
cohort[, dose_quartile := cut(daily_dose_main,
  breaks = quantile(daily_dose_main, c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE),
  labels = c("Q1","Q2","Q3","Q4"),
  include.lowest = TRUE)]

# Route IV (indicateur de VPA IV = souvent ICU)
cohort[, vpa_iv := as.integer(!is.na(has_iv) & has_iv)]

# =============================================================================
# 4. Score de sévérité approché
# =============================================================================
log_msg("4. Score de sévérité approché...")

# Score simplifié (somme pondérée de marqueurs de sévérité)
cohort[, severity_score := (
  as.integer(in_icu) * 3 +
  as.integer(!is.na(sepsis) & sepsis) * 2 +
  as.integer(crea_aki) +
  as.integer(alb_baseline_low) +
  as.integer(emergency) +
  as.integer(liver_baseline_elevated)
)]

# =============================================================================
# 5. Variables de modification d'effet prédéfinies
# =============================================================================
log_msg("5. Effect modifiers prédéfinis...")

# Déjà présents : comed_topiramate, comed_carbapenem, etc.
# Albumine basse (< 3)
cohort[, em_low_albumin := alb_baseline_low]
# Insuffisance rénale
cohort[, em_aki := crea_aki]
# Âge > 65
cohort[, em_elderly := as.integer(age > 65)]
# ICU
cohort[, em_icu := in_icu]
# Plaquettes proche seuil
cohort[, em_plt_borderline := as.integer(!is.na(plt_baseline) & plt_baseline < 150)]

# =============================================================================
# 6. Gestion des données manquantes
# =============================================================================
log_msg("6. Analyse des données manquantes...")

# Variables à imputer
vars_to_impute <- c("alb_baseline", "crea_baseline", "alt_baseline", "ast_baseline",
                     "bili_baseline", "weight_kg_baseline", "nh3_baseline")

miss_summary <- missing_summary(cohort[, vars_to_impute, with=FALSE])
log_msg("Résumé missing data:")
print(miss_summary)
safe_write_csv(miss_summary, file.path(TABLES_DIR, "missing_data_summary.csv"), overwrite=TRUE)

# Imputation multiple avec MICE pour les variables avec < 40% manquant
vars_mice <- vars_to_impute[sapply(vars_to_impute, function(v) {
  pct <- mean(is.na(cohort[[v]])) * 100
  pct < 40
})]

log_msg(sprintf("Variables pour MICE: %s", paste(vars_mice, collapse=", ")))

# Préparation données pour MICE — prédicteurs simples pour éviter la singularité
mice_predictors <- c("age", "female", "in_icu", "emergency", "high_dose")
mice_data_cols <- unique(c(mice_predictors, vars_mice))
mice_data_cols <- mice_data_cols[mice_data_cols %in% names(cohort)]
mice_data <- cohort[, mice_data_cols, with=FALSE]

# Convertir les logiques en numériques
for (col in names(mice_data)) {
  if (is.logical(mice_data[[col]])) mice_data[, (col) := as.integer(get(col))]
}

set.seed(RANDOM_SEED)
imp <- tryCatch(
  mice(as.data.frame(mice_data), m=5, method="pmm",
       maxit=10, seed=RANDOM_SEED, printFlag=FALSE),
  error = function(e) {
    log_msg(sprintf("MICE erreur, utilisation norm: %s", e$message), level="WARN")
    # Essayer avec méthode norm
    mice(as.data.frame(mice_data), m=5, method="norm",
         maxit=5, seed=RANDOM_SEED, printFlag=FALSE)
  }
)
log_msg("Imputation MICE complète")

# Extraire les données imputées (dataset complet = 1) et marquer les originaux manquants
imp_complete <- as.data.table(complete(imp, action=1))

# Remplacer les colonnes imputées dans la cohorte
for (v in vars_mice) {
  if (v %in% names(imp_complete)) {
    imp_col <- imp_complete[[v]]
    imp_name <- paste0(v, "_imp")
    miss_name <- paste0(v, "_missing")
    cohort[, (miss_name) := as.integer(is.na(get(v)))]
    cohort[, (imp_name) := imp_col]
  } else {
    # Imputation simple : médiane pour les variables non dans MICE
    med_val <- median(cohort[[v]], na.rm=TRUE)
    miss_name <- paste0(v, "_missing")
    imp_name  <- paste0(v, "_imp")
    cohort[, (miss_name) := as.integer(is.na(get(v)))]
    cohort[[imp_name] ] <- ifelse(is.na(cohort[[v]]), med_val, cohort[[v]])
  }
}

# Sauvegarder les 5 datasets imputés
for (i in 1:5) {
  arrow::write_parquet(as.data.frame(complete(imp, action=i)),
                        file.path(DATA_INTERMEDIATE, sprintf("mice_imp_%d.parquet", i)))
}
saveRDS(imp, file.path(DATA_INTERMEDIATE, "mice_object.rds"))
log_msg("Objets MICE sauvegardés")

# =============================================================================
# 7. Sélection des covariables finales pour le modèle de propension
# =============================================================================
log_msg("7. Jeu de covariables final...")

# Variables de propension
propensity_vars <- c(
  # Démographie
  "age", "female", "race_simple",
  # Sévérité
  "in_icu", "sepsis", "emergency", "pre_hosp_days_winsorized",
  "epilepsy_dx", "bipolar_dx",
  # Biologie baseline
  "plt_baseline", "alb_baseline", "crea_baseline", "alt_baseline", "ast_baseline",
  "liver_baseline_elevated",
  # Co-médications
  "comed_topiramate", "comed_carbapenem", "comed_heparin", "comed_linezolid",
  "comed_steroid", "comed_antipsy", "comed_other_aed",
  # VPA
  "vpa_iv"
)

# Vérifier quelles variables sont disponibles
propensity_vars_avail <- propensity_vars[propensity_vars %in% names(cohort)]
propensity_vars_missing <- propensity_vars[!propensity_vars %in% names(cohort)]
if (length(propensity_vars_missing) > 0) {
  log_msg(sprintf("Variables non disponibles dans la cohorte: %s",
                  paste(propensity_vars_missing, collapse=", ")), level="WARN")
}

log_msg(sprintf("Variables de propension disponibles: %d", length(propensity_vars_avail)))

# Sauvegarder la liste des variables
saveRDS(propensity_vars_avail, file.path(DATA_INTERMEDIATE, "propensity_vars.rds"))

# =============================================================================
# 8. Sauvegarde finale
# =============================================================================
log_msg("8. Sauvegarde features...")
safe_write_parquet(cohort, file.path(DATA_FINAL, "features_cohort.parquet"), overwrite=TRUE)
safe_write_csv(cohort, file.path(DATA_FINAL, "features_cohort.csv"), overwrite=TRUE)

timer_end(t0, "Étape 6 totale")

# Rapport
feature_report <- sprintf(
"# Feature Engineering — Étude VPA/MIMIC-IV
**Date :** %s

## Variables créées

### Variables d'exposition
- `daily_dose_main` : dose journalière max dans les 72h (mg/j)
- `high_dose` : binaire (> 1000 mg/j) — **variable de traitement principale**
- `dose_cat_750/1000/1250/1500` : seuils alternatifs pour sensibilité
- `dose_per_kg` : dose rapportée au poids (mg/kg/j)
- `log_dose` : log(dose + 1) pour analyses continues
- `vpa_iv` : administration IV vs PO

### Variables biologiques baselines
- `alb_baseline_low` : albumine < 3 g/dL
- `crea_aki` : créatinine >= 1.5 mg/dL
- `liver_baseline_elevated` : ALT ou AST > 40 UI/L
- `plt_baseline_100_150` : plaquettes proches du seuil

### Effect modifiers prédéfinis
- `em_low_albumin`, `em_aki`, `em_elderly`, `em_icu`, `em_plt_borderline`

### Score de sévérité approché
- `severity_score` : somme pondérée (ICU×3, sepsis×2, AKI×1, alb_low×1, urgence×1)

## Imputation des données manquantes

Variables soumises à MICE (5 imputations) :
%s

5 datasets imputés sauvegardés : `data_intermediate/mice_imp_[1-5].parquet`

## Couverture des variables de propension

Variables disponibles : %d / %d

## Anti-leakage temporel

TOUTES les variables de covariables sont définies AVANT le time zero :
- Biologie baseline : fenêtre [-7j, 0j]
- Co-médications : prescriptions actives avant time zero
- Variables fixes : âge, sexe, race
- Contexte : ICU avant time zero, diagnostic à l'admission

Aucune variable post-time-zero n'est incluse dans le jeu de covariables.
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(sprintf("- `%s` : %.0f%% manquant", vars_mice,
                sapply(vars_mice, function(v) mean(is.na(cohort[[v]]))*100)),
        collapse="\n"),
  length(propensity_vars_avail),
  length(propensity_vars)
)

writeLines(feature_report, file.path(REPORTS_DIR, "06_feature_engineering.md"))
log_msg("=== ÉTAPE 6 TERMINÉE ===")
