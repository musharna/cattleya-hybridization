#!/usr/bin/env Rscript
# What predicts breeders' use of a species (pre-registration: prereg/U_species_use.md).
#
# What: prereg C's test-3 model (Poisson T1, exposure offset, phylogenetic + observation random effects)
# with WCVP predictors added: native to Brazil, wet tropical climate, lithophyte, standardized log range size.
# Reports each predictor's rate ratio (posterior median, 95% CrI). Tree, data, sampler and fitting come
# from signal_corrected.R, sourced unchanged (its modes do not match "use"/"usesim", so it only sets up).
# Usage, from analysis/ (C_TREE/C_TAG/C_ITER/C_THIN/C_ADAPT as for signal_corrected.R):
#   Rscript species_use.R usesim 5 planted SHARD   # planted rate ratio 2 on every predictor, lambda 0.3
#   Rscript species_use.R usesim 5 null SHARD      # no predictor effect, lambda 0.3
#   Rscript species_use.R use                      # real counts
# Output: ../data/species_use_{real,sim}<tag>*.json

source("signal_corrected.R")
stopifnot(mode %in% c("use", "usesim"))
cv <- read.csv("../data/species_covariates.csv", stringsAsFactors = FALSE)
cv <- cv[match(d$Binomial, cv$Binomial), ]
stopifnot(identical(cv$Binomial, d$Binomial), !anyNA(cv))
X <- data.frame(brazil = cv$brazil, wet = cv$wet, litho = cv$litho, z_areas = as.numeric(scale(log(cv$n_areas))))
preds <- names(X)
F_U <- T1 ~ transferred + brazil + wet + litho + z_areas + offset(logE) + (1 | gr(species, cov = A)) + (1 | obs)

rr <- function(fit) {
  dr <- as_draws_df(fit)
  out <- lapply(c("transferred", preds), function(p) {
    v <- exp(dr[[paste0("b_", p)]])
    c(rr = median(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)))
  })
  c(setNames(out, c("transferred", preds)), list(h2 = h2(fit)))
}

# sim_counts() from signal_corrected.R with a fixed-effect shift added on the log scale (same draws, same order)
sim_u <- function(l, shift) {
  y <- log1p(d$primary_T1)
  bm <- fastBM(tree, sig2 = 1)[tree$tip.label]; bm <- (bm - mean(bm)) / sd(bm)
  z <- sqrt(l) * bm + sqrt(1 - l) * rnorm(nrow(d))
  mu <- (d$logE - mean(d$logE)) + mean(y) + sd(y) * z + shift
  rpois(nrow(d), exp(mu))
}

base <- data.frame(logE = d$logE, transferred = d$transferred, X)
if (mode == "usesim") {
  eff <- args[3]; shard <- as.integer(args[4])
  stopifnot(eff %in% c("planted", "null"), shard >= 1)
  set.seed(1000 * shard + if (eff == "planted") 50 else 60)
  shift <- if (eff == "planted") as.vector(scale(as.matrix(X), scale = FALSE) %*% rep(log(2), length(preds))) else 0
  dat <- transform(base, T1 = sim_u(0.3, shift))
  f <- brms_fit(dat, F_U)
  res <- list(rr(f))
  for (i in seq_len(nsim - 1)) {
    dat$T1 <- sim_u(0.3, shift)
    res[[i + 1]] <- rr(brms_fit(dat, F_U, prev = f))
  }
  write_json(list(nsim = nsim, effect = eff, shard = shard, reps = res),
             sprintf("../data/species_use_sim%s_%s_s%d.json", sfx, eff, shard), digits = 6, auto_unbox = TRUE)
}
if (mode == "use") {
  f <- brms_fit(transform(base, T1 = d$primary_T1), F_U)
  out <- rr(f)
  print(out)
  write_json(out, sprintf("../data/species_use_real%s.json", sfx), digits = 6, auto_unbox = TRUE)
}
