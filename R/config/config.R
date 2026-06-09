# =============================================================================
# config.R
# Global configuration: paths, constants, parameters
# =============================================================================

library(here)

# Répertoire racine du projet
PROJECT_ROOT <- here::here()

# Chemins MIMIC-IV
MIMIC_ROOT  <- file.path(dirname(PROJECT_ROOT), "mimic-iv-3.1")
MIMIC_HOSP  <- file.path(MIMIC_ROOT, "hosp")
MIMIC_ICU   <- file.path(MIMIC_ROOT, "icu")

# Chemins projet
DATA_RAW          <- file.path(PROJECT_ROOT, "data_raw")
DATA_INTERMEDIATE <- file.path(PROJECT_ROOT, "data_intermediate")
DATA_FINAL        <- file.path(PROJECT_ROOT, "data_final")
SQL_DIR           <- file.path(PROJECT_ROOT, "sql")
REPORTS_DIR       <- file.path(PROJECT_ROOT, "reports")
OUTPUTS_DIR       <- file.path(PROJECT_ROOT, "outputs")
TABLES_DIR        <- file.path(OUTPUTS_DIR, "tables")
FIGURES_DIR       <- file.path(OUTPUTS_DIR, "figures")
DIAGNOSTICS_DIR   <- file.path(OUTPUTS_DIR, "diagnostics")
MODELS_DIR        <- file.path(OUTPUTS_DIR, "models")
LOGS_DIR          <- file.path(OUTPUTS_DIR, "logs")

# =============================================================================
# Paramètres de l'étude
# =============================================================================

# Fenêtres temporelles
BASELINE_WINDOW_DAYS    <- 7   # jours avant time zero pour baseline labs
MAX_FOLLOWUP_DAYS       <- 30  # jours de suivi post-initiation
GRACE_PERIOD_DAYS       <- 2   # tolérance pour assignment

# Définitions des outcomes (thrombopénie)
THROMBOCYTOPENIA_ABS    <- 100   # G/L — seuil absolu
THROMBOCYTOPENIA_REL30  <- 0.30  # baisse relative 30%
THROMBOCYTOPENIA_REL50  <- 0.50  # baisse relative 50%

# Hyperammoniémie
AMMONIA_THRESHOLD_HIGH  <- 55   # µmol/L (valeur standard adulte)
AMMONIA_THRESHOLD_MOD   <- 35   # seuil modéré

# Hépatotoxicité
ALT_THRESHOLD_RATIO     <- 3    # x ULN (upper limit of normal)
AST_THRESHOLD_RATIO     <- 3
ALT_ABSOLUTE_HIGH       <- 120  # UI/L approximatif 3xULN moyen
AST_ABSOLUTE_HIGH       <- 120

# Exposition valproate
VPA_HIGH_DOSE_THRESHOLD_MG  <- 1000  # mg/24h
VPA_HIGH_LEVEL_THRESHOLD_UG <- 75    # µg/mL

# Itération random seed
RANDOM_SEED <- 20250430

# =============================================================================
# Identifiants labevents clés (à confirmer par feasibility)
# =============================================================================
# Ces codes seront vérifiés dans 01_feasibility.R
LAB_CODES <- list(
  # Plaquettes
  platelets = c(51265),          # Platelet Count
  # Ammoniaque
  ammonia   = c(50971),          # Ammonia
  # Transaminases
  alt       = c(50861),          # ALT (SGPT)
  ast       = c(50878),          # AST (SGOT)
  bilirubin_total = c(50885),    # Bilirubin Total
  albumin   = c(50862),          # Albumin
  creatinine = c(50912),         # Creatinine
  # Valproate (niveau sérique)
  valproate_level = c(50819, 52742)  # Valproic Acid — à confirmer
)

# =============================================================================
# Poids dans chartevents ICU
# =============================================================================
WEIGHT_ITEM_IDS <- c(226512, 226531, 224639)  # Admission Weight, Daily Weight

# Logging
LOG_TIMESTAMP_FORMAT <- "%Y-%m-%d %H:%M:%S"

cat("Config chargée. PROJECT_ROOT:", PROJECT_ROOT, "\n")
cat("MIMIC_ROOT:", MIMIC_ROOT, "\n")
