# =============================================================================
# utils.R
# Utility functions: logging, I/O, SMD, missing data
# =============================================================================

library(data.table)
library(lubridate)

# Logging horodaté
log_msg <- function(..., level = "INFO") {
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", ts, "] [", level, "] ", paste(..., sep = " "))
  cat(msg, "\n")
  invisible(msg)
}

# Sauvegarde sécurisée (ne remplace pas sans backup)
safe_write_parquet <- function(dt, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    backup <- sub("\\.parquet$", paste0("_backup_", format(Sys.time(), "%Y%m%d%H%M%S"), ".parquet"), path)
    file.copy(path, backup)
    log_msg("Backup créé:", backup, level = "WARN")
  }
  arrow::write_parquet(as.data.frame(dt), path)
  log_msg("Fichier écrit:", path)
}

safe_write_csv <- function(dt, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    backup <- sub("\\.csv$", paste0("_backup_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv"), path)
    file.copy(path, backup)
    log_msg("Backup créé:", backup, level = "WARN")
  }
  data.table::fwrite(as.data.table(dt), path)
  log_msg("CSV écrit:", path)
}

# Lecture parquet
read_parquet_dt <- function(path) {
  log_msg("Lecture:", path)
  as.data.table(arrow::read_parquet(path))
}

# Résumé rapide d'un data.table
dt_summary <- function(dt, name = "") {
  cat("\n===", name, "===\n")
  cat("Lignes:", nrow(dt), " | Colonnes:", ncol(dt), "\n")
  cat("Colonnes:", paste(names(dt), collapse = ", "), "\n\n")
}

# Vérification que les répertoires existent
ensure_dirs <- function(...) {
  dirs <- c(...)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      log_msg("Répertoire créé:", d)
    }
  }
}

# Calcul de SMD (Standardized Mean Difference)
smd <- function(x, group, na.rm = TRUE) {
  x1 <- x[group == 1]
  x0 <- x[group == 0]
  if (na.rm) {
    x1 <- x1[!is.na(x1)]
    x0 <- x0[!is.na(x0)]
  }
  mu1 <- mean(x1)
  mu0 <- mean(x0)
  s1  <- var(x1)
  s0  <- var(x0)
  pooled_sd <- sqrt((s1 + s0) / 2)
  if (pooled_sd == 0) return(NA_real_)
  (mu1 - mu0) / pooled_sd
}

# Résumé de missing data
missing_summary <- function(dt) {
  miss <- sapply(dt, function(x) sum(is.na(x)))
  pct  <- miss / nrow(dt) * 100
  data.table(
    variable    = names(miss),
    n_missing   = miss,
    pct_missing = round(pct, 2)
  )[order(-n_missing)]
}

# Winsorisation d'une variable numérique
winsorize <- function(x, lower = 0.01, upper = 0.99) {
  q <- quantile(x, c(lower, upper), na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

# Timer simple
timer_start <- function() proc.time()
timer_end   <- function(start, label = "") {
  elapsed <- proc.time() - start
  log_msg(sprintf("%s | Durée: %.1f sec", label, elapsed["elapsed"]))
}
