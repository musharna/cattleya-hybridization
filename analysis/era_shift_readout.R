#!/usr/bin/env Rscript
# Readout for prereg V: planted-truth gate, then (if the real fit exists) H2_shift, the mini-catt mechanism
# check and one figure of per-species shifts in tree order. Rules are prereg V's, applied mechanically.
# Usage, from analysis/: Rscript era_shift_readout.R
# Output: ../data/era_shift_readout.json, ../figures_regenerated/deeper/V_era_shift.png
suppressPackageStartupMessages({library(jsonlite); library(ape)})
source("house_theme.R")
rd <- function(g) do.call(rbind, lapply(Sys.glob(g), function(f) fromJSON(f)$reps))
hi_ <- rd("../data/era_shift_sim_mol_l0.6_s*.json"); lo_ <- rd("../data/era_shift_sim_mol_l0_s*.json")
stopifnot(nrow(hi_) == 10, nrow(lo_) == 10)
cnt <- function(m) { c_ <- m[m[, 4] < 1.01, , drop = FALSE]; c(pass = sum(c_[, 2] >= 0.05), n = nrow(c_)) }
p6 <- cnt(hi_); p0 <- cnt(lo_)
invalid <- p0[["pass"]] > 2; underpowered <- p6[["pass"]] < 8
cat(sprintf("planted lambda_s 0.6: %d of %d converged pass; lambda_s 0: %d of %d; invalid %s, underpowered %s\n",
            p6[["pass"]], p6[["n"]], p0[["pass"]], p0[["n"]], invalid, underpowered))
out <- list(planted = p6, null = p0, invalid = invalid, underpowered = underpowered)

real_f <- "../data/era_shift_real_mol.json"
if (file.exists(real_f)) {
  r <- fromJSON(real_f)
  h <- unlist(r$h2_shift); sh <- unlist(r$shift)
  founders <- paste0("Cattleya_", c("coccinea", "walkeriana", "aclandiae", "intermedia"))
  stopifnot(all(founders %in% names(sh)))
  mech <- sum(sh[founders] > 0)
  verdict <- if (h[4] >= 1.01) "unconverged" else if (mech < 3) "not interpretable (mechanism check failed)" else
    if (h[2] >= 0.05 && !invalid && !underpowered) "supported" else
    sprintf("not supported (%s)", paste(c(if (h[2] < 0.05) "lower bound < 0.05", if (invalid) "invalid",
                                            if (underpowered) "underpowered"), collapse = ", "))
  cat(sprintf("H2_shift %.3f [%.4f, %.3f], Rhat %.4f; era effect (late) %.2f; mini-catt founders with positive shift %d of 4 (%s); verdict: %s\n",
              h[1], h[2], h[3], h[4], r$b_late[1], mech, paste(sprintf("%.2f", sh[founders]), collapse = ", "), verdict))
  out$real <- list(h2_shift = h, b_late = r$b_late, founders = sh[founders], mech = mech, verdict = verdict)
  tree <- read.tree("../data/molecular_combined.nwk"); tree <- ladderize(tree)
  ord <- tree$tip.label[tree$edge[tree$edge[, 2] <= length(tree$tip.label), 2]]
  ord <- ord[ord %in% names(sh)]
  df <- data.frame(sp = factor(sub("Cattleya_", "C. ", ord), levels = sub("Cattleya_", "C. ", ord)), shift = sh[ord],
                   founder = ord %in% founders)
  p <- ggplot(df, aes(shift, sp, colour = founder)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = pal_null) + geom_point(size = 1.8) +
    scale_colour_manual(values = c(`FALSE` = pal_unknown, `TRUE` = pal_obs), labels = c("other", "mini-catt founder (A)"), name = NULL) +
    labs(x = "Change in primary-hybrid rate after 1960, relative to the average species (log scale)", y = NULL,
         title = "Which species gained or lost breeders' use after 1960",
         subtitle = sprintf("Species in tree order. Phylogenetic share of the change: %.2f (%.3f-%.2f)", h[1], h[2], h[3])) +
    theme_house(10) + theme(axis.text.y = element_text(size = 5.5), legend.position = "bottom")
  dir.create("../figures_regenerated/deeper", recursive = TRUE, showWarnings = FALSE)
  save_house(p, "../figures_regenerated/deeper/V_era_shift.png", w = 7, h = 11)
}
write_json(out, "../data/era_shift_readout.json", digits = 6, auto_unbox = TRUE)
