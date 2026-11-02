#!/usr/bin/env Rscript
# Readout for prereg U: planted-truth gates per predictor, then (if the real fit exists) the label per
# predictor and one figure of rate ratios. Rules are prereg U's, applied mechanically.
# Usage, from analysis/: Rscript species_use_readout.R
# Output: ../data/species_use_readout.json, ../figures_regenerated/deeper/U_species_use.png
suppressPackageStartupMessages(library(jsonlite))
source("house_theme.R")
preds <- c("brazil", "wet", "litho", "z_areas")
rd <- function(f) fromJSON(f, simplifyVector = FALSE)$reps
pl <- unlist(lapply(Sys.glob("../data/species_use_sim_mol_planted_s*.json"), rd), recursive = FALSE)
nu <- unlist(lapply(Sys.glob("../data/species_use_sim_mol_null_s*.json"), rd), recursive = FALSE)
stopifnot(length(pl) == 10, length(nu) == 10)
conv <- function(reps) Filter(function(r) r$h2[[4]] < 1.01, reps)
plc <- conv(pl); nuc <- conv(nu)
gate <- do.call(rbind, lapply(preds, function(p) {
  rec <- sum(sapply(plc, function(r) r[[p]][[2]] > 1))
  flag <- sum(sapply(nuc, function(r) r[[p]][[2]] > 1 || r[[p]][[3]] < 1))
  data.frame(predictor = p, recovered = rec, n_planted = length(plc), flagged = flag, n_null = length(nuc),
             underpowered = rec < 8, invalid = flag > 2)
}))
print(gate)
out <- list(gate = gate, unconverged_planted = length(pl) - length(plc), unconverged_null = length(nu) - length(nuc))

real_f <- "../data/species_use_real_mol.json"
if (file.exists(real_f)) {
  r <- fromJSON(real_f, simplifyVector = FALSE)
  out$real_converged <- r$h2[[4]] < 1.01
  if (!out$real_converged) message("REAL FIT UNCONVERGED (Rhat ", r$h2[[4]], "): labels below are reported as unconverged")
  tab <- do.call(rbind, lapply(c(preds, "transferred"), function(p) {
    v <- unlist(r[[p]]); g <- gate[gate$predictor == p, ]
    lab <- if (p == "transferred") "as in C" else
      if (g$invalid) "invalid" else
      if (v[2] > 1 || v[3] < 1) "matters" else
      if (v[2] >= 0.5 && v[3] <= 2 && !g$underpowered) "no strong effect" else "unresolved"
    data.frame(predictor = p, rr = v[1], lo = v[2], hi = v[3], label = lab)
  }))
  print(tab); cat(sprintf("H2 with predictors %.3f [%.4f, %.3f], Rhat %.4f\n", r$h2[[1]], r$h2[[2]], r$h2[[3]], r$h2[[4]]))
  out$real <- tab; out$h2 <- unlist(r$h2)
  nice <- c(brazil = "native to Brazil", wet = "wet tropical climate", litho = "lithophyte",
            z_areas = "range size (per SD of log areas)", transferred = "moved into Cattleya")
  tab$what <- factor(sprintf("%s: %s", nice[tab$predictor], tab$label), levels = rev(sprintf("%s: %s", nice[tab$predictor], tab$label)))
  p <- ggplot(tab, aes(rr, what)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = pal_null) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", colour = pal_cat[1]) +
    geom_point(size = 3, colour = pal_cat[1]) + scale_x_log10(breaks = c(0.5, 0.7, 1, 1.4, 2)) +
    labs(x = "Rate ratio, primary hybrids per year (median, 95% CrI)", y = NULL,
         title = "What predicts how often breeders used a species",
         subtitle = "91 species; allows for time since description and the phylogeny") +
    theme_house() + theme(plot.title.position = "plot")
  dir.create("../figures_regenerated/deeper", recursive = TRUE, showWarnings = FALSE)
  save_house(p, "../figures_regenerated/deeper/U_species_use.png", w = 9, h = 3.8)
}
write_json(out, "../data/species_use_readout.json", digits = 6, auto_unbox = TRUE)
