#!/usr/bin/env Rscript
# Moved species before vs after the 2008 genus moves (pre-registration: prereg/T_moved_period.md).
#
# What: Poisson counts of primaries per species x period, offset log(exposure), moved x post
# interaction, phylogenetic + species + observation random effects on the molecular tree.
# ratio = exp(moved:post) = rate ratio after / before. sim mode plants a known ratio first.
# Usage, from analysis/:
#   Rscript moved_period.R sim constant|naming NREP SHARD   # planted truth (gate)
#   Rscript moved_period.R real                              # only after the gate passes
# Output: ../data/moved_period_<mode>[_<scenario>_s<shard>].json

suppressPackageStartupMessages({ library(ape); library(phytools); library(jsonlite); library(brms) })
args <- commandArgs(TRUE)
mode <- args[1]
tree <- read.tree("../data/molecular_combined.nwk")
pc <- read.csv("../data/species_period_counts.csv", stringsAsFactors = FALSE)
tr <- read.csv("../data/species_traits.csv", stringsAsFactors = FALSE)
pc$species <- gsub(" ", "_", pc$Binomial)
pc <- pc[match(tree$tip.label, pc$species), ]
stopifnot(!anyNA(pc$species), is.ultrametric(tree, tol = 1e-6))
A <- vcv(tree, corr = TRUE)
sd_y <- sd(log1p(tr$primary_T1[match(pc$Binomial, tr$Binomial)]))

long <- function(n_pre, n_post, e_pre, e_post) {
  d <- rbind(data.frame(species = pc$species, moved = pc$moved, post = 0, n = n_pre, E = e_pre),
             data.frame(species = pc$species, moved = pc$moved, post = 1, n = n_post, E = e_post))
  d <- d[d$E >= 1, ]
  d$logE <- log(d$E); d$species_iid <- d$species; d$obs <- seq_len(nrow(d))
  d
}
e_pre <- pc$cutoff - pc$descr_year
e_post <- 2025 - pc$cutoff
F_T <- n ~ moved * post + offset(logE) + (1 | gr(species, cov = A)) + (1 | species_iid) + (1 | obs)
fit <- function(d, form = F_T, prev = NULL) {
  if (is.null(prev))
    brm(form, data = d, data2 = list(A = A), family = poisson(), chains = 4, iter = 4000, cores = 4,
        seed = 1, refresh = 0, silent = 2, control = list(adapt_delta = 0.95))
  else update(prev, newdata = d, data2 = list(A = A), seed = 1, refresh = 0, silent = 2)
}
readout <- function(f, term = "moved:post") {
  s <- fixef(f)[term, c("Estimate", "Q2.5", "Q97.5")]
  r <- exp(unname(s))
  c(ratio = r[1], lo = r[2], hi = r[3], rhat_max = max(rhat(f), na.rm = TRUE))
}
verdict <- function(r, placebo_ok = TRUE) {
  if (r[["rhat_max"]] >= 1.01) return("unconverged")
  if (r[["hi"]] < 1 && placebo_ok) return("naming")
  if (r[["lo"]] <= 1 && r[["hi"]] >= 1 && r[["lo"]] >= 0.5 && r[["hi"]] <= 2) return("appeal")
  "cant_tell"
}

if (mode == "sim") {
  scen <- args[2]; nsim <- as.integer(args[3]); shard <- as.integer(args[4])
  stopifnot(scen %in% c("constant", "naming"))
  set.seed(1000 * shard + if (scen == "naming") 7 else 3)
  tot <- c(sum(pc$n_pre), sum(pc$n_post))  # pooled totals only
  b_post <- if (scen == "naming") 0 else log(2.5)
  sim1 <- function() {
    bm <- fastBM(tree, sig2 = 1)[tree$tip.label]; bm <- (bm - mean(bm)) / sd(bm)
    z <- sd_y * (sqrt(0.3) * bm + sqrt(0.7) * rnorm(nrow(pc)))
    n <- list()
    for (p in 1:2) {
      E <- pmax(if (p == 1) e_pre else e_post, 1e-9)
      rest <- log(E) + pc$moved * (if (p == 1) log(2.5) else b_post) + z + rnorm(nrow(pc), 0, 0.3)
      a <- log(tot[p] / sum(exp(rest)))
      n[[p]] <- rpois(nrow(pc), exp(a + rest))
    }
    long(n[[1]], n[[2]], e_pre, e_post)
  }
  f0 <- fit(sim1())
  res <- list(readout(f0))
  for (i in seq_len(nsim - 1)) res[[i + 1]] <- readout(fit(sim1(), prev = f0))
  res <- as.data.frame(do.call(rbind, res))
  res$verdict <- apply(res, 1, function(r) verdict(as.list(r)))
  res$scenario <- scen
  print(res)
  write_json(list(scenario = scen, nsim = nsim, reps = res),
             sprintf("../data/moved_period_sim_%s_s%d.json", scen, shard), digits = 6)
}

if (mode == "real") {
  d <- long(pc$n_pre, pc$n_post, e_pre, e_post)
  f <- fit(d)
  dp <- long(pc$n_p1985, pc$n_p2008, 1985 - pc$descr_year, rep(2008 - 1985, nrow(pc)))
  fp <- fit(dp, prev = f)
  d14 <- long(pc$n_pre14, pc$n_post14, 2014 - pc$descr_year, rep(2025 - 2014, nrow(pc)))
  f14 <- fit(d14, prev = f)
  ds <- d; ds$sat <- ds$post * log1p(pc$n_pre[match(ds$species, pc$species)])
  fs <- fit(ds, update(F_T, . ~ . + sat))
  pl <- readout(fp)
  main <- readout(f)
  out <- list(rows = nrow(d), rows_dropped = 2 * nrow(pc) - nrow(d), main = main, placebo_1985 = pl,
              cutoff_2014 = readout(f14), saturation = readout(fs),
              moved_pre = exp(fixef(f)["moved", c("Estimate", "Q2.5", "Q97.5")]),
              verdict = verdict(as.list(main), placebo_ok = pl[["lo"]] <= 1 && pl[["hi"]] >= 1))
  print(out)
  write_json(out, "../data/moved_period_real.json", digits = 6, auto_unbox = TRUE)
}
