#!/usr/bin/env Rscript
# Readout for prereg W: planted-truth gates (power, false positives, CI coverage), then, if the real run
# exists, the frozen label, the Python second-route check, and one figure (observed / expected primaries by
# distance decile, with the permutation null band). Rules are prereg W's, applied mechanically.
# Usage, from analysis/: Rscript species_pairs_readout.R
# Output: ../data/species_pairs_readout.json, ../figures_regenerated/deeper/W_species_pairs.png
suppressPackageStartupMessages(library(jsonlite))
rd <- function(rr, sdv) {
  fs <- Sys.glob(sprintf("../data/species_pairs_sim_rr%s_sd%s_*.json", rr, sdv))
  m <- do.call(rbind, lapply(fs, function(f) fromJSON(f)$reps)); colnames(m) <- c("beta", "se", "rr", "lo", "hi", "p", "n_species"); m
}
arm <- function(rr, sdv, n_expected) {
  m <- rd(rr, sdv); stopifnot(nrow(m) == n_expected)
  hit <- if (rr < 1) m[, "p"] < 0.05 & m[, "beta"] < 0 else if (rr > 1) m[, "p"] < 0.05 & m[, "beta"] > 0 else m[, "p"] < 0.05
  c(rr = rr, sd = sdv, n = nrow(m), detected = sum(hit), covered = sum(m[, "lo"] <= rr & m[, "hi"] >= rr),
    median_rr = median(m[, "rr"]))
}
g <- rbind(arm(0.67, 1, 100), arm(1.5, 1, 100), arm(1, 1, 100))
sens <- rbind(arm(1, 0.5, 50), arm(1, 1.5, 50), arm(0.67, 0.5, 50), arm(0.67, 1.5, 50))
print(g); print(sens)
up_close <- g[1, "detected"] < 80; up_far <- g[2, "detected"] < 80
invalid <- g[3, "detected"] > 8
uncal <- any(g[, "covered"] < 85)
cat(sprintf("gates: underpowered close %s, distant %s; invalid %s; CI uncalibrated %s\n", up_close, up_far, invalid, uncal))
out <- list(gates = g, sensitivity = sens, underpowered_close = up_close, underpowered_distant = up_far,
            invalid = invalid, ci_uncalibrated = uncal)

real_f <- "../data/species_pairs_real.json"
if (file.exists(real_f)) {
  r <- fromJSON(real_f); m <- unlist(r$main); stopifnot(!is.null(names(m)))
  lab <- if (invalid) "invalid" else if (m[["p"]] < 0.05 && m[["rr"]] < 1) "closer pairs crossed more" else
    if (m[["p"]] < 0.05 && m[["rr"]] > 1) "more distant pairs crossed more" else
    if (m[["lo"]] >= 0.67 && m[["hi"]] <= 1.5 && !up_close && !up_far) "no strong preference" else "unresolved"
  py <- fromJSON("../data/species_pairs_check_py.json")
  stopifnot(abs(py$beta - m[["beta"]]) < 0.01)
  cat(sprintf("real: RR %.3f per SD (95%% CI %.3f-%.3f), permutation p %.3f, %d species; label: %s; Python beta %.5f vs R %.5f\n",
              m[["rr"]], m[["lo"]], m[["hi"]], m[["p"]], m[["n_species"]], lab, py$beta, m[["beta"]]))
  out$real <- list(main = m, label = lab, python_beta = py$beta, no_natural = r$no_natural, its = r$its,
                   boot_range = apply(r$boot, 2, range), era = r$era)

  # figure: observed / expected primaries by distance decile; expected from the species-only model (beta = 0);
  # band = the same ratio under 199 label permutations
  commandArgs <- function(...) "lib"; source("species_pairs.R"); rm(commandArgs); source("house_theme.R")
  set.seed(1)
  f0 <- fit_fast(tomat(n_obs), Dm * 0)                       # species effects only
  kp <- f0$keep[pr[, 1]] & f0$keep[pr[, 2]]
  mu0 <- f0$mu; nk <- n_obs[kp]; dk <- dv[kp]
  dec <- function(d) cut(d, quantile(dv[kp], 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  oe <- function(d) { b <- dec(d); tapply(nk, b, sum) / tapply(mu0, b, sum) }
  # mu0 and nk are in upper-triangle order of the kept species; map pair order to it
  ks <- which(f0$keep); up <- which(upper.tri(matrix(0, length(ks), length(ks))), arr.ind = TRUE)
  ord <- match(paste(ks[up[, 1]], ks[up[, 2]]), paste(pr[kp, 1], pr[kp, 2])); stopifnot(!anyNA(ord))
  nk <- n_obs[kp][ord]; dk <- dv[kp][ord]
  obs <- oe(dk)
  nul <- replicate(199, { o <- sample(S); oe(Dm[o, o][cbind(pr[kp, 1], pr[kp, 2])][ord]) })
  df <- data.frame(decile = 1:10, obs = obs, lo = apply(nul, 1, quantile, 0.025), hi = apply(nul, 1, quantile, 0.975))
  p <- ggplot(df, aes(decile)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = pal_null, alpha = 0.35) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = pal_null) +
    geom_line(aes(y = obs), colour = pal_cat[1]) + geom_point(aes(y = obs), colour = pal_cat[1], size = 2.5) +
    scale_x_continuous(breaks = 1:10, labels = c("closest", 2:9, "most distant")) +
    labs(x = "Species pairs by phylogenetic distance (deciles)", y = "Observed / expected primaries",
         title = "Are closely related species crossed more often?",
         subtitle = sprintf("Expected allows for each species' popularity and how long both were known.\nGrey: 95%% range when distances are shuffled. Rate ratio per SD of distance %.2f (%.2f-%.2f), p %.3f",
                            m[["rr"]], m[["lo"]], m[["hi"]], m[["p"]])) +
    theme_house() + theme(plot.title.position = "plot")
  dir.create("../figures_regenerated/deeper", recursive = TRUE, showWarnings = FALSE)
  save_house(p, "../figures_regenerated/deeper/W_species_pairs.png", w = 9, h = 4.5)
  out$figure <- df
}
write_json(out, "../data/species_pairs_readout.json", digits = 6, auto_unbox = TRUE)
