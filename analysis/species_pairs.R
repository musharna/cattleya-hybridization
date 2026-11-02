#!/usr/bin/env Rscript
# Are registered primary crosses more or less common between closely related species? (prereg/W_species_pairs.md)
#
# What: Poisson pseudo-likelihood gravity model over all 4,095 pairs of the 91 molecular-tree species,
# n_ij ~ beta * d_ij + a_i + a_j + offset(log E_ij). beta per SD of cophenetic distance; 95% CI from a
# dyadic-robust sandwich; p from permuting species labels on the distance matrix (QAP).
# Usage, from analysis/:
#   Rscript species_pairs.R selftest                 # margin identity, known-answer recovery, timing
#   Rscript species_pairs.R sim RR SD START N NPERM  # planted truth, RR in {1, 0.67, 1.5}
#   Rscript species_pairs.R real                     # real counts + sensitivities + era interaction
# Output: ../data/species_pairs_sim_rr<RR>_sd<SD>_<START>.json, ../data/species_pairs_real.json
suppressPackageStartupMessages({library(ape); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "selftest"

tr <- read.tree("../data/molecular_combined.nwk")
sp <- read.csv("../data/species_traits.csv", stringsAsFactors = FALSE)
sp <- sp[match(gsub("_", " ", tr$tip.label), sp$Binomial), ]
stopifnot(nrow(sp) == 91, !anyNA(sp$Binomial))
S <- nrow(sp); tips <- tr$tip.label
pr <- t(combn(S, 2)); K <- nrow(pr)                            # pair k = species pr[k,1], pr[k,2]
X <- matrix(0, K, S); X[cbind(1:K, pr[, 1])] <- 1; X[cbind(1:K, pr[, 2])] <- 1
oy <- sp$orig_year
logE <- log(pmax(1, 2025 - pmax(oy[pr[, 1]], oy[pr[, 2]])))

dist_vec <- function(tree) {
  D <- cophenetic(tree)[tips, tips]
  v <- D[pr]; (v - mean(v)) / sd(v)
}
dv <- dist_vec(tr)

pc <- read.csv("../data/species_pair_counts.csv", stringsAsFactors = FALSE)
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")
kk <- key(sp$Binomial[pr[, 1]], sp$Binomial[pr[, 2]])
pc <- pc[key(pc$sp1, pc$sp2) %in% kk, ]
n_obs <- numeric(K); n_obs[match(key(pc$sp1, pc$sp2), kk)] <- pc$n
stopifnot(sum(n_obs) == 1157)  # printed by derive_pair_counts.py before any modelling

# Reference engine: glm.fit with species effects as two-member indicators (column 1 dropped), plus distance
fit_glm <- function(n, d, off = logE, Xs = X) {
  M <- cbind(d = d, Xs[, -1, drop = FALSE], int = 1)
  f <- glm.fit(M, n, family = poisson(), offset = off, control = list(epsilon = 1e-10, maxit = 100))
  stopifnot(f$converged)
  unname(f$coefficients[1])
}
# Fast engine, same MLE: species effects by damped Jacobi proportional fitting on the S x S matrix (each
# species' fitted total = observed total), beta by Newton on the distance score; alternate to convergence.
# Species with zero observed total have no finite MLE effect and carry no information about beta: dropped.
tomat <- function(v) { m <- matrix(0, S, S); m[pr] <- v; m + t(m) }
Em <- exp(tomat(logE)); diag(Em) <- 0
fit_fast <- function(nm, Dmat, keep = rowSums(nm) > 0, tol = 1e-9) {
  nm <- nm[keep, keep]; Dk <- Dmat[keep, keep]; Ek <- Em[keep, keep]; N <- rowSums(nm)
  a <- log(N / rowSums(Ek)) / 2; b <- 0; up <- upper.tri(nm)
  for (it in 1:20000) {
    W <- Ek * exp(b * Dk)
    for (j in 1:5) a <- (a + log(N / (W %*% exp(a))[, 1])) / 2
    mu <- exp(outer(a, a, "+")) * W
    g <- sum(((nm - mu) * Dk)[up]); h <- sum((mu * Dk^2)[up])
    b <- b + if (h > 0) g / h else 0  # h = 0 only for a distance-free (species-only) fit
    if (abs(g) < tol * sum(N) && max(abs(rowSums(mu) - N)) < tol * sum(N)) break
  }
  stopifnot(it < 20000)
  list(beta = b, mu = mu[up], n = nm[up], keep = keep, it = it)
}
# dyadic-robust variance of beta (pairs sharing a species are correlated): meat = sum_s S_s S_s' - sum_k u_k u_k'
se_dyadic <- function(f, d) {
  ks <- which(f$keep); Xk <- matrix(0, length(f$n), length(ks))
  pk <- which(upper.tri(matrix(0, length(ks), length(ks))), arr.ind = TRUE)
  Xk[cbind(seq_along(f$n), pk[, 1])] <- 1; Xk[cbind(seq_along(f$n), pk[, 2])] <- 1
  M <- cbind(d, Xk[, -1, drop = FALSE], 1)
  u <- M * (f$n - f$mu)
  B <- solve(crossprod(M * sqrt(f$mu)))
  meat <- crossprod(crossprod(Xk, u)) - crossprod(u)
  sqrt((B %*% meat %*% B)[1, 1])
}
Dstd <- function(tree) { D <- cophenetic(tree)[tips, tips]; v <- D[pr]; (D - mean(v)) / sd(v) }
Dm <- Dstd(tr)
summ <- function(n, nperm, Dmat = Dm) {
  nm <- tomat(n); f <- fit_fast(nm, Dmat)
  dk <- Dmat[f$keep, f$keep][upper.tri(nm[f$keep, f$keep])]
  se <- se_dyadic(f, dk)
  b <- replicate(nperm, { o <- sample(S); fit_fast(nm, Dmat[o, o], keep = f$keep)$beta })
  c(beta = f$beta, se = se, rr = exp(f$beta), lo = exp(f$beta - 1.96 * se), hi = exp(f$beta + 1.96 * se),
    p = (1 + sum(abs(b) >= abs(f$beta))) / (1 + nperm), n_species = sum(f$keep))
}

# planted truth: species effects from the published per-species totals only; c scales the expected total to 1157
a_sp <- log((sp$primary_T1 + 0.5) / pmax(1, 2025 - oy)); a_sp <- a_sp - mean(a_sp)
sim_n <- function(rr, sdev) {
  eta <- logE + a_sp[pr[, 1]] + a_sp[pr[, 2]] + log(rr) * dv
  c0 <- log(1157 / sum(exp(eta + sdev^2 / 2)))
  rpois(K, exp(c0 + eta + rnorm(K, 0, sdev)))
}

if (mode == "selftest") {
  set.seed(1)
  n1 <- sim_n(0.5, 0); f1 <- fit_fast(tomat(n1), Dm)
  cat(sprintf("known answer (RR 0.5, no noise): fast RR %.4f in %d iterations, total %d\n", exp(f1$beta), f1$it, sum(n1)))
  stopifnot(abs(exp(f1$beta) - 0.5) < 0.1)
  # the two engines must agree on the same data (glm.fit on the kept pairs only)
  for (r in 1:3) {
    n2 <- sim_n(c(1, 0.67, 1.5)[r], 1); f2 <- fit_fast(tomat(n2), Dm)
    kp <- f2$keep[pr[, 1]] & f2$keep[pr[, 2]]
    bg <- fit_glm(n2[kp], dv[kp], logE[kp], X[kp, f2$keep])
    cat(sprintf("engine check %d: fast %.6f, glm %.6f, species kept %d\n", r, f2$beta, bg, sum(f2$keep)))
    stopifnot(abs(f2$beta - bg) < 1e-5)
  }
  t0 <- proc.time()[["elapsed"]]; s1 <- summ(sim_n(1, 1), 199); dt <- proc.time()[["elapsed"]] - t0
  print(round(s1, 4)); cat(sprintf("one summ() with 199 permutations: %.1f s\n", dt))
}
if (mode == "sim") {
  rr <- as.numeric(args[2]); sdev <- as.numeric(args[3]); st <- as.integer(args[4]); nn <- as.integer(args[5]); np <- as.integer(args[6])
  stopifnot(rr %in% c(1, 0.67, 1.5), sdev %in% c(0.5, 1, 1.5))
  base <- c(`1` = 3000, `0.67` = 4000, `1.5` = 5000)[[as.character(rr)]]
  res <- t(sapply(st:(st + nn - 1), function(r) {
    set.seed(base + r + if (sdev == 1) 0 else round(sdev * 10) * 1000)
    out <- summ(sim_n(rr, sdev), np); cat(r, round(out, 4), "\n"); out
  }))
  write_json(list(rr = rr, sd = sdev, start = st, nperm = np, reps = res),
             sprintf("../data/species_pairs_sim_rr%s_sd%s_%d.json", rr, sdev, st), digits = 6)
}
