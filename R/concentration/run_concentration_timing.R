# =============================================================================
# run_concentration_timing.R
# Conc-6 — Trough timing reclassification analysis
# =============================================================================

# Run from the project root directory (where R/ and analysis.qmd are located).
source("R/config/config.R")
source("R/utils/utils.R")

library(data.table); library(arrow); library(ggplot2); library(patchwork)

log_msg("=== STARTING ===")
t0 <- timer_start()
ensure_dirs(TABLES_DIR, FIGURES_DIR)

# =============================================================================
# 1. Chargement données
# =============================================================================
vpa_lev <- as.data.table(read_parquet(file.path(DATA_INTERMEDIATE, "vpa_levels_wk1.parquet")))
vpa_pres <- as.data.table(read_parquet(file.path(DATA_INTERMEDIATE, "vpa_prescriptions_raw.parquet")))

vpa_lev[, charttime := as.POSIXct(charttime)]
vpa_pres[, starttime := as.POSIXct(starttime)]
vpa_pres[, stoptime := as.POSIXct(stoptime)]

cat("N mesures VPA levels:", nrow(vpa_lev), "\n")
cat("N prescriptions:", nrow(vpa_pres), "\n")

# =============================================================================
# 2. Construction du dataset avec heure + métadonnées prescription
# =============================================================================
vpa_lev[, hour_of_day := as.integer(format(charttime, "%H"))]

# Enrichir avec route et fréquence de la prescription contemporaine
setkey(vpa_lev, subject_id)
setkey(vpa_pres, subject_id)

enrich <- rbindlist(lapply(1:nrow(vpa_lev), function(i) {
  sid <- vpa_lev$subject_id[i]
  ct  <- vpa_lev$charttime[i]
  pres <- vpa_pres[subject_id == sid & starttime <= ct & stoptime >= ct]
  if (nrow(pres) == 0) pres <- vpa_pres[subject_id == sid][which.min(abs(as.numeric(difftime(starttime, ct, units="hours"))))]
  if (nrow(pres) == 0) return(data.table(route="Inconnu", doses_per_24h=NA_real_, interval_h=NA_real_))
  data.table(
    route = pres$route[1],
    doses_per_24h = pres$doses_per_24_hrs[1],
    interval_h = ifelse(!is.na(pres$doses_per_24_hrs[1]) & pres$doses_per_24_hrs[1] > 0,
                        24/pres$doses_per_24_hrs[1], NA_real_)
  )
}))

dt <- cbind(vpa_lev[, .(subject_id, charttime, hour_of_day, value_ugml)], enrich)
dt <- dt[!is.na(value_ugml)]
cat("N mesures avec valeur:", nrow(dt), "\n")

# =============================================================================
# 3. Quatre classifications alternatives
# =============================================================================

# Classification 1 — STRICTE (classification actuelle V2)
dt[, class_strict := fcase(
  hour_of_day %in% 6:9,           "Trough probable (6h-9h)",
  hour_of_day %in% c(12,13,20,21),"Post-dose probable",
  default="Timing indéterminé"
)]

# Classification 2 — ELARGIE (fenêtre matinale 4h-11h comme trough plausible)
# Rationale : administration VPA souvent à 22h-23h; résiduel attendu à 4h-11h
dt[, class_elargie := fcase(
  hour_of_day %in% 4:11,          "Trough plausible (4h-11h)",
  hour_of_day %in% c(12,13,20,21),"Post-dose plausible",
  default="Timing indéterminé"
)]

# Classification 3 — PHARMACOLOGIQUE (croise route + fréquence)
# IV : perfusion continue = steady-state (aucune classification trough/peak valide)
# PO/OG BID (2/j) : trough attendu ~12h après dernière prise
# PO/OG QD (1/j)  : trough attendu ~20-24h après dernière prise
# PO/OG TID (3/j) : intervalle ~8h
dt[, class_pharma := fcase(
  route %in% c("IV"),                                "Perfusion IV (steady-state attendu)",
  doses_per_24h %in% c(1,2) & hour_of_day %in% 4:11,"Trough BID/QD plausible (4h-11h)",
  doses_per_24h %in% c(3,4) & hour_of_day %in% 5:9, "Trough TID/QID plausible (5h-9h)",
  hour_of_day %in% c(12,13,20,21),                  "Post-dose plausible",
  default="Timing indéterminé"
)]

# Classification 4 — TRES CONSERVATIVE (heure de prélèvement très précise)
dt[, class_conservative := fcase(
  hour_of_day %in% 6:8,           "Trough probable (6h-8h seulement)",
  hour_of_day %in% c(12,13,20,21),"Post-dose probable",
  default="Timing indéterminé"
)]

# =============================================================================
# 4. Tableau de comparaison des 4 classifications
# =============================================================================
classify_summary <- function(col_name, label) {
  col <- dt[[col_name]]
  n_indet <- sum(grepl("indéterminé", col, ignore.case=TRUE))
  n_trough <- sum(grepl("rough|steady", col, ignore.case=TRUE))
  n_postdose <- sum(grepl("ost-dose|peak", col, ignore.case=TRUE))
  data.table(
    classification = label,
    n_total = nrow(dt),
    n_trough_plausible = n_trough,
    pct_trough = round(n_trough/nrow(dt)*100, 1),
    n_postdose = n_postdose,
    pct_postdose = round(n_postdose/nrow(dt)*100, 1),
    n_indeterminate = n_indet,
    pct_indeterminate = round(n_indet/nrow(dt)*100, 1)
  )
}

reclassification_summary <- rbindlist(list(
  classify_summary("class_strict",      "1. Stricte (6h-9h = trough) — V2"),
  classify_summary("class_elargie",     "2. Elargie (4h-11h = trough plausible)"),
  classify_summary("class_pharma",      "3. Pharmacologique (route × fréquence)"),
  classify_summary("class_conservative","4. Conservative (6h-8h uniquement)")
))
log_msg("Table de reclassification:")
print(reclassification_summary)

# Changements de classification : combien passent de "indéterminé" (strict) à "trough plausible" (élargi)?
n_newly_trough <- sum(grepl("indéterminé", dt$class_strict) &
                      grepl("Trough|steady", dt$class_elargie, ignore.case=TRUE))
n_total_indet_strict <- sum(grepl("indéterminé", dt$class_strict))
pct_reclassified <- round(n_newly_trough / nrow(dt) * 100, 1)
cat(sprintf("\nReclassification élargie : %d/%d (%.1f%%) passent de 'indéterminé' à 'trough plausible'\n",
            n_newly_trough, n_total_indet_strict, pct_reclassified))

# Stratification par route
route_timing <- dt[, .(
  n = .N,
  pct_trough_strict = round(mean(grepl("Trough prob", class_strict))*100, 1),
  pct_trough_elargi = round(mean(grepl("Trough plausible", class_elargie))*100, 1),
  pct_indet_strict  = round(mean(grepl("indéterminé", class_strict, ignore.case=TRUE))*100, 1),
  pct_indet_elargi  = round(mean(grepl("indéterminé", class_elargie, ignore.case=TRUE))*100, 1)
), by=route][order(-n)]
log_msg("Timing par route:"); print(route_timing)

# Sauvegarder la table de reclassification
safe_write_csv(reclassification_summary,
               file.path(TABLES_DIR, "concentration_timing_reclassification.csv"), overwrite=TRUE)

# =============================================================================
# 5. Figure enrichie — histogramme horaire avec 3 classifications + zones
# =============================================================================

# Préparation longue pour facets
dt_long <- rbindlist(list(
  dt[, .(hour_of_day, classification_type = "1. Stricte (6h-9h)",
         timing_class = class_strict, value_ugml, route)],
  dt[, .(hour_of_day, classification_type = "2. Elargie (4h-11h)",
         timing_class = class_elargie, value_ugml, route)],
  dt[, .(hour_of_day, classification_type = "3. Pharmacologique (route×fréq.)",
         timing_class = class_pharma, value_ugml, route)]
))

# Palette commune
pal <- c(
  "Trough probable (6h-9h)"            = "#2196F3",
  "Trough plausible (4h-11h)"          = "#1976D2",
  "Trough BID/QD plausible (4h-11h)"   = "#1565C0",
  "Trough TID/QID plausible (5h-9h)"   = "#0D47A1",
  "Trough probable (6h-8h seulement)"  = "#42A5F5",
  "Perfusion IV (steady-state attendu)" = "#78909C",
  "Post-dose probable"                  = "#F44336",
  "Post-dose plausible"                 = "#E53935",
  "Timing indéterminé"                  = "#9E9E9E"
)

p_main <- ggplot(dt_long, aes(x=hour_of_day, fill=timing_class)) +
  geom_histogram(bins=24, color="white", linewidth=0.15) +
  scale_fill_manual(values=pal, name="Classification") +
  scale_x_continuous(breaks=seq(0,23,3), labels=sprintf("%02dh",seq(0,23,3))) +
  facet_wrap(~classification_type, ncol=1, scales="free_y") +
  labs(
    title="Distribution horaire des prélèvements VPA — Analyse de sensibilité des classifications",
    subtitle=sprintf("N=%d mesures J0-J7 | Cluster matinal 5h-7h compatible avec prélèvements résiduels", nrow(dt)),
    x="Heure du prélèvement", y="Nombre de mesures"
  ) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=9),
        strip.background=element_rect(fill="gray95"))

# Graphique des proportions par classification
prop_data <- reclassification_summary[, .(
  classification,
  Trough=pct_trough, Indéterminé=pct_indeterminate, `Post-dose`=pct_postdose
)]
prop_long <- melt(prop_data, id.vars="classification",
                  variable.name="Catégorie", value.name="Pourcentage")
prop_long[, classification := factor(classification, levels=rev(reclassification_summary$classification))]

p_prop <- ggplot(prop_long, aes(x=classification, y=Pourcentage, fill=Catégorie)) +
  geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%.0f%%", Pourcentage)), position=position_stack(vjust=0.5),
            size=3, color="white", fontface="bold") +
  scale_fill_manual(values=c("Trough"="#2196F3","Indéterminé"="#9E9E9E","Post-dose"="#F44336")) +
  coord_flip() +
  labs(title="Proportion par classification",
       x=NULL, y="% des prélèvements") +
  theme_bw(base_size=11) + theme(legend.position="bottom")

p_combined <- p_main / p_prop + plot_layout(heights=c(3,1))
ggsave(file.path(FIGURES_DIR, "vpa_concentration_timing_analysis_v2.png"),
       p_combined, width=12, height=13, dpi=150)
log_msg("Figure timing v2 sauvegardée")

# =============================================================================
# 6. Texte de conclusion nuancée
# =============================================================================
pct_indet_strict <- reclassification_summary[grepl("Stricte", classification), pct_indeterminate]
pct_indet_elargi <- reclassification_summary[grepl("Elargie", classification), pct_indeterminate]
pct_trough_elargi <- reclassification_summary[grepl("Elargie", classification), pct_trough]

cat(sprintf("\n=== CONCLUSION TIMING V3 ===\n"))
cat(sprintf("Classification stricte (V2): indéterminé = %.1f%%\n", pct_indet_strict))
cat(sprintf("Classification élargie (V3): indéterminé = %.1f%% (trough plausible = %.1f%%)\n",
            pct_indet_elargi, pct_trough_elargi))
cat(sprintf("Reclassifiés de indéterminé → trough plausible: %d mesures (%.1f%% du total)\n",
            n_newly_trough, pct_reclassified))

conclusion_text <- sprintf(
  "La classification initiale (stricte, 6h-9h) produisait %.1f%% de timing indéterminé. ",
  pct_indet_strict
)
conclusion_text <- paste0(conclusion_text,
  sprintf("Une classification élargie (4h-11h) compatible avec les pratiques réelles de prélèvement en ICU reclasse %.1f%% des prélèvements vers 'trough plausible', ", pct_reclassified),
  sprintf("réduisant le timing indéterminé à %.1f%%.", pct_indet_elargi)
)

if (pct_indet_elargi > 40) {
  conclusion_text <- paste0(conclusion_text,
    " Même avec cette classification plus souple, le timing reste indéterminé pour une fraction substantielle, soutenant le maintien de la qualification 'exploratoire' pour l'analyse concentration-based — mais sans en invalider totalement les conclusions directionnelles.")
} else {
  conclusion_text <- paste0(conclusion_text,
    " Sous la classification pharmacologique élargie, la majorité des prélèvements est compatible avec des concentrations résiduelles, assouplissant partiellement la critique initiale.")
}
cat("\n", conclusion_text, "\n")

# Sauvegarder conclusion
writeLines(c(
  "# Conclusion — Reclassification Timing VPA V3",
  paste0("Date: ", Sys.Date()),
  "",
  "## Résultats des 4 classifications alternatives",
  "",
  capture.output(print(reclassification_summary)),
  "",
  sprintf("## Reclassification clé : %d mesures (%.1f%%) passent de 'indéterminé' (classification stricte) à 'trough plausible' (classification élargie 4h-11h)",
          n_newly_trough, pct_reclassified),
  "",
  "## Conclusion",
  conclusion_text,
  "",
  "## Position recommandée dans le rapport",
  "- Analyse concentration-based demeure EXPLORATOIRE",
  "- Mais la qualification 'timing 58.3% indéterminé' mérite nuance : fondée sur une règle conservatrice",
  "- Avec une fenêtre réaliste (4h-11h), le timing indéterminé se réduit sensiblement",
  "- Ce constat ne valide pas une inférence causale robuste (ESS=3.2% reste critique)",
  "- Formulation recommandée : voir appendix_concentration_timing_reappraisal.md"
), file.path(TABLES_DIR, "timing_reappraisal_conclusion.txt"))

log_msg(sprintf("=== TIMING REAPPRAISAL COMPLET (%.1fs) ===", timer_stop(t0)))
