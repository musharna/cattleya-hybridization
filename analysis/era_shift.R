#!/usr/bin/env Rscript
# Is the pre/post-1960 change in species' primary-hybrid rates clustered on the phylogeny?
# (pre-registration: prereg/V_era_shift.md)
#
# What: one row per species x era (early <= 1959, late >= 1960), Poisson with exposure offset, an era effect,
# `transferred`, and phylogenetic + independent random intercepts and `late` slopes. H2_shift = phylogenetic
# share of the slope variance. Tree, data, sampler settings and h2-style summaries follow signal_corrected.R,
# sourced unchanged (its modes do not match "era"/"erasim", so it only sets up).
# Usage, from analysis/ (C_TREE/C_TAG/C_ITER/C_THIN/C_ADAPT as for signal_corrected.R):
#   Rscript era_shift.R erasim 5 0.6 SHARD   # planted phylogenetic slope structure lambda_s 0.6
#   Rscript era_shift.R erasim 5 0 SHARD     # none
#   Rscript era_shift.R era                  # real counts
# Output: ../data/era_shift_{real,sim}<tag>*.json

source("signal_corrected.R")
stopifnot(mode %in% c("era", "erasim"))
ec <- read.csv("../data/species_era_counts.csv", stringsAsFactors = FALSE)
ec <- ec[match(d$Binomial, ec$Binomial), ]
stopifnot(identical(ec$Binomial, d$Binomial))
early <- !is.na(ec$n_early)
rows <- rbind(
  data.frame(i = which(early), late = 0L, cnt = ec$n_early[early], logE = log(ec$E_early[early])),
  data.frame(i = seq_len(nrow(d)), late = 1L, cnt = ec$n_late, logE = log(ec$E_late)))
rows$species <- tree$tip.label[rows$i]; rows$sp <- rows$species
rows$transferred <- d$transferred[rows$i]; rows$obs <- seq_len(nrow(rows))
F_V <- cnt ~ late + transferred + offset(logE) + (1 + late || gr(species, cov = A)) + (1 + late || sp) + (1 | obs)

fit_v <- function(dat, prev = NULL) {
  suppressPackageStartupMessages(library(brms))
  A <- vcv(tree, corr = TRUE)
  if (is.null(prev))
    brm(F_V, data = dat, data2 = list(A = A), family = poisson(), chains = 4, iter = iter, thin = thin,
        cores = 4, seed = 1, refresh = 0, silent = 2, control = list(adapt_delta = adapt))
  else update(prev, newdata = dat, data2 = list(A = A), seed = 1, refresh = 0, silent = 2)
}
summ_v <- function(fit) {
  dr <- as_draws_df(fit)
  sp_ <- dr$sd_species__late^2; si_ <- dr$sd_sp__late^2
  h <- sp_ / (sp_ + si_)
  # per-species shift = deviation from the average change (Amendment 1): phylogenetic + independent late effects
  shift <- sapply(tree$tip.label, function(s) median(dr[[sprintf("r_species[%s,late]", s)]] + dr[[sprintf("r_sp[%s,late]", s)]]))
  stopifnot(length(sp_) > 0, length(si_) > 0, length(shift) == length(tree$tip.label), all(is.finite(shift)))
  list(h2_shift = c(median = median(h), lo = unname(quantile(h, 0.025)), hi = unname(quantile(h, 0.975)),
                    rhat_max = max(brms::rhat(fit), na.rm = TRUE)),
       b_late = unname(quantile(dr$b_late, c(0.5, 0.025, 0.975))), shift = as.list(shift))
}

# Simulation as in sim_counts() (intercepts at planted lambda 0.3), plus a standardized per-species late shift
# with phylogenetic share lambda_s, and the era effect calibrated from the real totals only
b_era <- log((sum(ec$n_late) / sum(ec$E_late)) / (sum(ec$n_early, na.rm = TRUE) / sum(ec$E_early, na.rm = TRUE)))
sim_v <- function(ls) {
  y <- log1p(d$primary_T1)
  bm <- fastBM(tree, sig2 = 1)[tree$tip.label]; bm <- (bm - mean(bm)) / sd(bm)
  z <- sqrt(0.3) * bm + sqrt(0.7) * rnorm(nrow(d))
  b2 <- fastBM(tree, sig2 = 1)[tree$tip.label]; b2 <- (b2 - mean(b2)) / sd(b2)
  s <- sqrt(ls) * b2 + sqrt(1 - ls) * rnorm(nrow(d)); s <- (s - mean(s)) / sd(s)
  mu <- (rows$logE - mean(rows$logE)) + mean(y) + sd(y) * z[rows$i] + rows$late * (b_era + 1.0 * s[rows$i])
  rpois(nrow(rows), exp(mu))
}

if (mode == "erasim") {
  ls <- as.numeric(args[3]); shard <- as.integer(args[4])
  stopifnot(ls %in% c(0, 0.6), shard >= 1)
  set.seed(1000 * shard + if (ls > 0) 70 else 80)
  dat <- transform(rows, cnt = sim_v(ls))
  f <- fit_v(dat)
  res <- list(summ_v(f)$h2_shift)
  for (k in seq_len(nsim - 1)) {
    dat$cnt <- sim_v(ls)
    res[[k + 1]] <- summ_v(fit_v(dat, prev = f))$h2_shift
  }
  write_json(list(nsim = nsim, lambda_s = ls, shard = shard, reps = do.call(rbind, res)),
             sprintf("../data/era_shift_sim%s_l%s_s%d.json", sfx, ls, shard), digits = 6)
}
if (mode == "era") {
  f <- fit_v(rows)  # fit first: summ_v() needs brms, which fit_v() loads (a lazy argument would not)
  out <- summ_v(f)
  print(out$h2_shift)
  write_json(out, sprintf("../data/era_shift_real%s.json", sfx), digits = 6, auto_unbox = TRUE)
}
