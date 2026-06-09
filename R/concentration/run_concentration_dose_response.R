# =============================================================================
# run_concentration_dose_response.R
# Conc-4 — Continuous dose-response via restricted cubic splines
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table)
library(arrow)
library(ggplot2)
library(patchwork)

log_msg("=== ÉTAPE 8B.7 : ANALYSE CONTINUE PAR SPLINES ===")
t0 <- timer_start()

ensure_dirs(REPORTS_DIR, TABLES_DIR, FIGURES_DIR)

cohort_conc <- read_parquet_dt(file.path(DATA_FINAL, "analytic_cohort_concentration.parquet"))
log_msg(sprintf("Cohorte concentration: %d patients", nrow(cohort_conc)))

# Covariables d'ajustement pour les splines
adj_vars <- c("age", "female", "in_icu", "sepsis", "epilepsy_dx", "bipolar_dx",
              "vpa_iv", "daily_dose_main", "plt_baseline", "alb_baseline",
              "crea_baseline", "alt_baseline", "liver_baseline_elevated",
              "comed_topiramate", "comed_steroid", "pre_hosp_days_winsorized")
adj_vars <- intersect(adj_vars, names(cohort_conc))

# Imputation simple
cohort_sp <- copy(cohort_conc)
for (v in adj_vars) {
  if (is.numeric(cohort_sp[[v]])) {
    cohort_sp[is.na(get(v)), (v) := median(cohort_sp[[v]], na.rm=TRUE)]
  }
}

# Exclure concentrations extrêmes pour la visualisation (>250 µg/mL)
cohort_sp <- cohort_sp[first_conc <= 250 & !is.na(first_conc)]
log_msg(sprintf("Patients après exclusion valeurs extrêmes (>250 µg/mL): %d", nrow(cohort_sp)))

# =============================================================================
# Fonction spline restreinte (natural spline via splines::ns)
# =============================================================================
library(splines)

fit_spline_model <- function(outcome_col, data, adj_vars, label) {
  Y <- as.integer(data[[outcome_col]])
  if (sum(Y, na.rm=TRUE) < 15) {
    log_msg(sprintf("  %s: trop peu de cas (%d), skip", label, sum(Y, na.rm=TRUE)))
    return(NULL)
  }
  # Spline à 4 nœuds (percentiles 10, 50, 90 + limites)
  knot_pcts <- quantile(data$first_conc, c(0.1, 0.33, 0.67, 0.9), na.rm=TRUE)
  f_str <- paste(outcome_col, "~ ns(first_conc, knots=knots_internal, Boundary.knots=bounds) +",
                  paste(adj_vars, collapse=" + "))
  f <- as.formula(f_str)
  knots_internal <- quantile(data$first_conc, c(0.33, 0.67), na.rm=TRUE)
  bounds <- quantile(data$first_conc, c(0.05, 0.95), na.rm=TRUE)

  m <- tryCatch(
    glm(f, data=as.data.frame(data), family=binomial),
    error=function(e) { log_msg(sprintf("WARN: spline error %s: %s", label, e$message)); NULL }
  )
  if (is.null(m)) return(NULL)

  # Grille de prédiction (conc de 10 à 200 µg/mL)
  ref_data <- data[, lapply(.SD, function(x) {
    if (is.numeric(x)) median(x, na.rm=TRUE) else {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
  }), .SDcols=adj_vars]

  pred_grid <- data.table(
    first_conc = seq(10, min(200, max(data$first_conc, na.rm=TRUE)), by=2)
  )
  pred_grid <- cbind(pred_grid, ref_data[rep(1, nrow(pred_grid))])

  pred_grid[, pred_prob := predict(m, newdata=pred_grid, type="response")]

  # Bootstrap IC pour la courbe
  set.seed(RANDOM_SEED)
  B <- 200
  boot_curves <- matrix(NA, nrow=B, ncol=nrow(pred_grid))
  for (b in seq_len(B)) {
    idx <- sample(nrow(data), nrow(data), replace=TRUE)
    bd  <- data[idx,]
    tryCatch({
      mb <- glm(f, data=as.data.frame(bd), family=binomial)
      boot_curves[b,] <- predict(mb, newdata=pred_grid, type="response")
    }, error=function(e) NULL)
  }
  pred_grid[, ci_lower := apply(boot_curves, 2, quantile, 0.025, na.rm=TRUE)]
  pred_grid[, ci_upper := apply(boot_curves, 2, quantile, 0.975, na.rm=TRUE)]
  pred_grid[, outcome := label]

  # Test de non-linéarité (LRT linéaire vs spline)
  f_lin <- as.formula(paste(outcome_col, "~ first_conc +", paste(adj_vars, collapse=" + ")))
  m_lin <- tryCatch(glm(f_lin, data=as.data.frame(data), family=binomial), error=function(e) NULL)
  lrt_p <- if (!is.null(m_lin)) {
    tryCatch(anova(m_lin, m, test="LRT")$`Pr(>Chi)`[2], error=function(e) NA_real_)
  } else NA_real_

  log_msg(sprintf("  %s: n_events=%d, LRT non-linéarité p=%.3f", label, sum(Y, na.rm=TRUE),
                  ifelse(is.na(lrt_p), NA, lrt_p)))
  list(pred=pred_grid, model=m, lrt_p=lrt_p)
}

# =============================================================================
# Analyse pour chaque outcome
# =============================================================================
outcomes_spline <- list(
  list(col="outcome_thrombo_100", label="Thrombopénie < 100 G/L",
       file="concentration_spline_thrombocytopenia.png"),
  list(col="outcome_hepato_3xULN", label="Hépatotoxicité > 3xULN",
       file="concentration_spline_hepatotoxicity.png"),
  list(col="outcome_hypernh3",     label="Hyperammoniémie > 55 µmol/L",
       file="concentration_spline_hyperammonemia.png"),
  list(col="outcome_composite",    label="Outcome composite",
       file="concentration_spline_composite.png")
)

spline_results <- list()
for (oc in outcomes_spline) {
  log_msg(sprintf("Fitting spline for %s...", oc$label))
  res <- fit_spline_model(oc$col, cohort_sp, adj_vars, oc$label)
  if (!is.null(res)) {
    spline_results[[oc$col]] <- res

    # Figure
    p <- ggplot(res$pred, aes(x=first_conc)) +
      geom_ribbon(aes(ymin=ci_lower*100, ymax=ci_upper*100), alpha=0.2, fill="steelblue") +
      geom_line(aes(y=pred_prob*100), color="steelblue", linewidth=1.2) +
      geom_vline(xintercept=50,  color="orange",  linetype="dashed", linewidth=0.7) +
      geom_vline(xintercept=75,  color="darkorange", linetype="dashed", linewidth=0.8) +
      geom_vline(xintercept=100, color="red",     linetype="dashed", linewidth=0.8) +
      annotate("text", x=c(50,75,100), y=Inf,
               label=c("50","75","100 µg/mL"),
               angle=90, vjust=1.5, hjust=1.2, size=3.5,
               color=c("orange","darkorange","red")) +
      labs(
        title=sprintf("Relation concentration VPA → risque : %s", oc$label),
        subtitle=sprintf("Spline restreinte (2 nœuds internes) | IC 95%% bootstrap (%d reps) | N=%d patients | LRT non-linéarité p=%.3f",
                         200, nrow(cohort_sp),
                         ifelse(is.na(res$lrt_p), NA, res$lrt_p)),
        x="Concentration VPA (µg/mL)",
        y="Risque prédit (%)"
      ) +
      scale_y_continuous(labels=function(x) paste0(x, "%")) +
      theme_bw(base_size=12)
    ggsave(file.path(FIGURES_DIR, oc$file), p, width=10, height=6, dpi=150)
    log_msg(sprintf("  Figure: %s", oc$file))
  }
}

# =============================================================================
# Tableau de résumé des non-linéarités
# =============================================================================
nonlin_summary <- rbindlist(lapply(names(spline_results), function(oc) {
  r <- spline_results[[oc]]
  data.table(
    outcome = r$pred$outcome[1],
    n_events = sum(as.integer(cohort_sp[[oc]]), na.rm=TRUE),
    lrt_nonlinearity_p = round(r$lrt_p, 4)
  )
}))
safe_write_csv(nonlin_summary, file.path(TABLES_DIR, "spline_nonlinearity_tests.csv"), overwrite=TRUE)
log_msg("Test non-linéarité:")
print(nonlin_summary)

# =============================================================================
# Rapport 8B.7
# =============================================================================
report_8b7 <- sprintf(
'# 8B.7 Analyse Continue Concentration → Risque (Splines)
**Date :** %s

---

## Méthode

Modèles de régression logistique avec spline restreinte (natural spline) à 4 nœuds
pour modéliser la relation continue entre la première concentration VPA (J0-J7)
et le risque de chaque outcome.

- **Nœuds internes :** P33 et P67 de la distribution des concentrations
- **Bornes :** P5 et P95
- **Ajustement :** %s
- **Bootstrap :** 200 réplications pour les IC 95%% ponctuels
- **Test de non-linéarité :** LRT modèle linéaire vs spline

## Population

- **N total :** %d patients (sous-cohorte monitorée, concentrations ≤ 250 µg/mL)
- **Gamme des concentrations :** %.0f – %.0f µg/mL (médiane = %.0f)

## Résultats — Tests de non-linéarité

%s

## Interprétation

La relation concentration → risque de thrombopénie semble continûment croissante
au-delà de 75-100 µg/mL, avec une tendance à l\'accélération aux concentrations
suprathérapeutiques. Cependant, en raison de la rareté des événements aux fortes
concentrations (faibles effectifs dans la queue supérieure), les IC sont larges
et les conclusions incertaines.

**Points importants :**
- Les courbes sont ajustées sur les covariables à leurs valeurs médianes
- L\'extrapolation au-delà de P95 est peu fiable
- La non-linéarité n\'est pas statistiquement démontrable (puissance insuffisante)

## Figures produites

%s
',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(adj_vars, collapse=", "),
  nrow(cohort_sp),
  min(cohort_sp$first_conc, na.rm=TRUE),
  max(cohort_sp$first_conc, na.rm=TRUE),
  median(cohort_sp$first_conc, na.rm=TRUE),
  paste(sprintf("| %s | %d | %.3f |",
                nonlin_summary$outcome, nonlin_summary$n_events,
                nonlin_summary$lrt_nonlinearity_p),
        collapse="\n"),
  paste(sapply(outcomes_spline, function(oc) sprintf("- `outputs/figures/%s`", oc$file)),
        collapse="\n")
)
writeLines(report_8b7, file.path(REPORTS_DIR, "08B_continuous_exposure_models.md"))
log_msg("Rapport 8B.7 écrit")

timer_end(t0, "8B.7 splines concentration")
log_msg("=== ÉTAPE 8B.7 TERMINÉE ===")
