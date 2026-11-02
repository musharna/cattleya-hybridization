#!/usr/bin/env Rscript
# Effort-controlled phylogenetic signal test (pre-registration: prereg/C_signal_corrected.md).
#
# What: rebuilds v1's species tree, reproduces v1 on `Descendants` (positive control), then tests
# primary-hybrid counts (T1) with Pagel's lambda, a description-era placebo, and a brms Poisson
# model with an exposure offset, a `transferred` effect and a phylogenetic random effect (H2).
# Planted-truth simulations on the same tree give the false-positive rate and power.
# Usage, from analysis/:
#   Rscript signal_corrected.R sim1 200     # test 1 on planted lambda 0 / 0.3 / 0.6
#   Rscript signal_corrected.R real         # control + tests 1-3 on the real traits
#   Rscript signal_corrected.R sim3 10      # test 3 on planted lambda 0 / 0.6
# Output: ../data/signal_corrected_<mode>.json, ../figures_regenerated/

suppressPackageStartupMessages({
  library(ape); library(phytools); library(jsonlite)
})
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "real"
nsim <- if (length(args) > 1) as.integer(args[2]) else 200
set.seed(1)
out_dir <- "../figures_regenerated"
dir.create(out_dir, showWarnings = FALSE)

# ---- tree: identical procedure to species_tree_signal.R --------------------
tree_lowestrank <- read.tree("../data/cattleya_subgeneric.nwk")
sp <- read.csv("../data/species_traits.csv", stringsAsFactors = FALSE)
v1 <- read.csv("../data/species_hybrid_counts.csv", stringsAsFactors = FALSE)
stopifnot(identical(sp$Binomial, v1$Binomial), identical(sp$LastSubgenus, v1$LastSubgenus),
          identical(as.numeric(sp$Descendants), as.numeric(v1$Descendants)))
sp$SpeciesName <- gsub(" ", "_", sp$Binomial)
rank_to_species <- split(sp$SpeciesName, sp$LastSubgenus)
tree_species <- tree_lowestrank
for (rank in names(rank_to_species)) {
  s <- rank_to_species[[rank]]
  idx <- which(tree_species$tip.label == rank)
  if (length(idx) == 0) next
  if (length(s) == 1) tree_species$tip.label[idx] <- s
  else tree_species <- bind.tree(tree_species, stree(length(s), tip.label = s, type = "star"), where = idx)
}
tree_species <- drop.tip(tree_species, intersect(tree_species$tip.label, tree_lowestrank$tip.label))
common <- intersect(tree_species$tip.label, sp$SpeciesName)
tree <- drop.tip(tree_species, setdiff(tree_species$tip.label, common))
tree$edge.length <- rep(1, nrow(tree$edge))
tree <- compute.brlen(tree, method = "Grafen")
# Amendment 2: C_TREE=<newick> replaces the v1 tree (ultrametric, tips Genus_species); C_TAG names outputs
tag <- Sys.getenv("C_TAG")
sfx <- if (nzchar(tag)) paste0("_", tag) else ""
if (nzchar(Sys.getenv("C_TREE"))) {
  tree <- read.tree(Sys.getenv("C_TREE"))
  stopifnot(nzchar(tag), is.ultrametric(tree, tol = 1e-6), all(tree$tip.label %in% sp$SpeciesName))
}
d <- sp[match(tree$tip.label, sp$SpeciesName), ]
d$E <- pmax(1, 2025 - d$orig_year)  # Amendment 5: original description, not the current name's year
d$logE <- log(d$E)
narrow <- c("Laelia", "Sophronitis", "Hoffmannseggella", "Sophronia")
d$transferred_narrow <- as.integer(d$basionym_genus %in% narrow)
cat(sprintf("tree tips %d; transferred %d (narrow %d)\n", nrow(d), sum(d$transferred), sum(d$transferred_narrow)))

lam <- function(x) { r <- phylosig(tree, setNames(x, tree$tip.label), method = "lambda", test = TRUE); c(lambda = r$lambda, p = r$P) }

sim_counts <- function(l) {
  y <- log1p(d$primary_T1)
  bm <- fastBM(tree, sig2 = 1)[tree$tip.label]; bm <- (bm - mean(bm)) / sd(bm)
  z <- sqrt(l) * bm + sqrt(1 - l) * rnorm(nrow(d))
  mu <- (d$logE - mean(d$logE)) + mean(y) + sd(y) * z
  rpois(nrow(d), exp(mu))
}

# Amendment 6: C_ITER/C_THIN/C_ADAPT override the sampler for the refits (defaults = Amendment 1 settings)
iter <- as.integer(Sys.getenv("C_ITER", "4000")); thin <- as.integer(Sys.getenv("C_THIN", "1"))
adapt <- as.numeric(Sys.getenv("C_ADAPT", "0.95"))
brms_fit <- function(dat, form, prev = NULL) {
  suppressPackageStartupMessages(library(brms))
  A <- vcv(tree, corr = TRUE)
  dat$species <- tree$tip.label; dat$obs <- seq_len(nrow(dat))
  if (is.null(prev))
    brm(form, data = dat, data2 = list(A = A), family = poisson(), chains = 4, iter = iter, thin = thin,
        cores = 4, seed = 1, refresh = 0, silent = 2, control = list(adapt_delta = adapt))
  else update(prev, newdata = dat, data2 = list(A = A), seed = 1, refresh = 0, silent = 2)
}
h2 <- function(fit) {
  dr <- as_draws_df(fit)
  s2p <- dr$sd_species__Intercept^2; s2o <- dr$sd_obs__Intercept^2
  h <- s2p / (s2p + s2o)
  c(median = median(h), lo = unname(quantile(h, 0.025)), hi = unname(quantile(h, 0.975)),
    rhat_max = max(brms::rhat(fit), na.rm = TRUE))
}
F_T <- T1 ~ transferred + offset(logE) + (1 | gr(species, cov = A)) + (1 | obs)

if (mode == "sim1") {
  res <- do.call(rbind, lapply(c(0, 0.3, 0.6), function(l) {
    t(sapply(seq_len(nsim), function(i) c(planted = l, lam(log1p(sim_counts(l))))))
  }))
  summ <- aggregate(cbind(lambda, sig = p < 0.05) ~ planted, data = as.data.frame(res), FUN = mean)
  print(summ)
  write_json(list(nsim = nsim, summary = summ, reps = as.data.frame(res)), sprintf("../data/signal_corrected_sim1%s.json", sfx), digits = 6)
}

if (mode == "real") {
  ctrl <- lam(d$Descendants)
  out <- list(
    control_v1_descendants = ctrl,
    control_pass = if (nzchar(sfx)) NA else unname(abs(ctrl["lambda"] - 0.337) <= 0.05 && ctrl["p"] < 1e-4),  # v1 tree only
    test1_T1 = lam(log1p(d$primary_T1)),
    test1_T1_K = { r <- phylosig(tree, setNames(log1p(d$primary_T1), tree$tip.label), method = "K", test = TRUE); c(K = r$K, p = r$P) },
    test1_T2 = lam(log1p(d$direct_T2)),
    placebo_logE = lam(d$logE))
  dat <- data.frame(T1 = d$primary_T1, logE = d$logE, transferred = d$transferred)
  f1 <- brms_fit(dat, F_T)
  f0 <- brms_fit(dat, T1 ~ offset(logE) + (1 | gr(species, cov = A)) + (1 | obs))
  dat_n <- transform(dat, transferred = d$transferred_narrow)
  fn <- brms_fit(dat_n, F_T, prev = f1)
  fx <- function(f) { s <- brms::fixef(f); s["transferred", c("Estimate", "Q2.5", "Q97.5")] }
  out$test3 <- list(with_transferred = h2(f1), without_transferred = h2(f0),
                    transferred_coef = fx(f1), narrow_with_transferred = h2(fn), narrow_coef = fx(fn))
  print(out)
  write_json(out, sprintf("../data/signal_corrected_real%s.json", sfx), digits = 6, auto_unbox = TRUE)
}

if (mode == "sim3") {
  # Rscript signal_corrected.R sim3 NREP LAMBDA SHARD  (shards run in parallel, distinct seeds)
  l <- as.numeric(args[3]); shard <- as.integer(args[4])
  set.seed(1000 * shard + round(10 * l))
  dat <- data.frame(T1 = sim_counts(l), logE = d$logE, transferred = d$transferred)
  f <- brms_fit(dat, F_T)
  res <- rbind(c(planted = l, h2(f)), t(sapply(seq_len(nsim - 1), function(i) {
    dat$T1 <- sim_counts(l); c(planted = l, h2(brms_fit(dat, F_T, prev = f)))
  })))
  print(res)
  write_json(list(nsim = nsim, reps = as.data.frame(res)),
             sprintf("../data/signal_corrected_sim3%s_l%s_s%d.json", sfx, l, shard), digits = 6)
}

if (mode == "boot") {
  # Rscript signal_corrected.R boot TREES.nwk FIRST LAST: test 3 on bootstrap trees FIRST..LAST, draws kept
  bt <- read.tree(args[2])[as.integer(args[3]):as.integer(args[4])]
  dat <- data.frame(T1 = d$primary_T1, logE = d$logE, transferred = d$transferred)
  hs <- list(); f <- NULL
  for (i in seq_along(bt)) {
    tree <- bt[[i]]
    stopifnot(setequal(tree$tip.label, d$SpeciesName))
    dd <- dat[match(tree$tip.label, d$SpeciesName), ]
    f <- brms_fit(dd, F_T, prev = f)
    dr <- as_draws_df(f)
    hs[[i]] <- list(h = dr$sd_species__Intercept^2 / (dr$sd_species__Intercept^2 + dr$sd_obs__Intercept^2),
                    rhat_max = max(brms::rhat(f), na.rm = TRUE))
  }
  write_json(list(trees = args[3]:args[4], h = lapply(hs, `[[`, "h"), rhat_max = sapply(hs, `[[`, "rhat_max")),
             sprintf("../data/signal_corrected_boot%s_%s_%s.json", sfx, args[3], args[4]), digits = 6)
}
