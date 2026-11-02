#!/usr/bin/env Rscript
# Real-data run for prereg W: main estimate (999 permutations), sensitivities (a) natural-hybrid pairs
# dropped, (b) 20 bootstrap trees, (c) ITS-only tree, and the secondary distance x era contrast.
# Loads species_pairs.R unchanged as a library (its mode switch sees "lib" and runs nothing).
# Usage, from analysis/: Rscript species_pairs_real.R
# Output: ../data/species_pairs_real.json, ../data/species_pairs_real_design.csv (for the Python second route)
commandArgs <- function(...) "lib"
source("species_pairs.R")
rm(commandArgs)
set.seed(20261029)
# W_DRY=1: plumbing test on simulated counts (planted RR 0.67) with 19 permutations; never touches real counts
dry <- Sys.getenv("W_DRY") == "1"; NP <- if (dry) 19 else 999
if (dry) {
  n_obs <- sim_n(0.67, 1); ne_ <- rbinom(K, n_obs, 0.5) * (pmax(oy[pr[, 1]], oy[pr[, 2]]) < 1960)
  pc <- data.frame(sp1 = sp$Binomial[pr[, 1]], sp2 = sp$Binomial[pr[, 2]], n = n_obs, n_early = ne_, n_late = n_obs - ne_)
  pc <- pc[pc$n > 0, ]
}
sfx <- if (dry) "_DRY" else ""

# fit_fast reads the exposure matrix Em from the global environment; with_E swaps it (0 = pair excluded)
with_E <- function(Enew, expr) { old <- Em; Em <<- Enew; on.exit(Em <<- old); expr }

out <- list(main = summ(n_obs, NP))
print(round(out$main, 4))

# second route: the kept-pair design for statsmodels
f0 <- fit_fast(tomat(n_obs), Dm); kp <- f0$keep[pr[, 1]] & f0$keep[pr[, 2]]
write.csv(data.frame(sp1 = sp$Binomial[pr[kp, 1]], sp2 = sp$Binomial[pr[kp, 2]], n = n_obs[kp], d = dv[kp], logE = logE[kp]),
          sprintf("../data/species_pairs_real_design%s.csv", sfx), row.names = FALSE)

# (a) pairs with a WCVP nothospecies excluded
np_ <- read.csv("../data/nothospecies_pairs.csv", stringsAsFactors = FALSE)
nat <- kk %in% key(np_$sp1, np_$sp2)
mask <- tomat(as.numeric(!nat)); n_a <- n_obs * !nat
out$no_natural <- with_E(Em * mask, summ(n_a, NP))
out$n_natural_pairs <- sum(nat); out$n_primaries_in_natural_pairs <- sum(n_obs[nat])

# (b) bootstrap trees, (c) ITS-only tree: point estimate and dyadic CI, no permutation
ci <- function(tree) {
  D <- Dstd(tree); f <- fit_fast(tomat(n_obs), D)
  dk <- D[f$keep, f$keep][upper.tri(D[f$keep, f$keep])]; se <- se_dyadic(f, dk)
  c(rr = exp(f$beta), lo = exp(f$beta - 1.96 * se), hi = exp(f$beta + 1.96 * se))
}
out$boot <- t(sapply(read.tree("../data/molecular_boot20.nwk"), ci))
out$its <- ci(read.tree("../data/molecular_ITS.nwk"))

# secondary: per-era fits (era-specific species effects and distance slope); contrast late - early,
# p by applying the same species permutation to both eras
pc2 <- pc[match(kk, key(pc$sp1, pc$sp2)), ]
ne <- ifelse(is.na(pc2$n_early), 0, pc2$n_early); nl <- ifelse(is.na(pc2$n_late), 0, pc2$n_late)
mo <- pmax(oy[pr[, 1]], oy[pr[, 2]])
E_e <- tomat(ifelse(mo < 1960, pmax(1, 1960 - pmax(mo, 1850)), 0))
E_l <- tomat(pmax(1, 2025 - pmax(mo, 1960)))
bad <- which(ne > 0 & mo >= 1960)  # an early grex of a pair not yet described before 1960 has no early exposure
if (length(bad)) stop("early primaries on pairs described >= 1960: ", paste(kk[bad], collapse = "; "))
fe <- with_E(E_e, fit_fast(tomat(ne), Dm)); fl <- with_E(E_l, fit_fast(tomat(nl), Dm))
dlt <- fl$beta - fe$beta
pb <- replicate(NP, { o <- sample(S)
  with_E(E_l, fit_fast(tomat(nl), Dm[o, o], keep = fl$keep)$beta) - with_E(E_e, fit_fast(tomat(ne), Dm[o, o], keep = fe$keep)$beta) })
out$era <- c(rr_early = exp(fe$beta), rr_late = exp(fl$beta), ratio_late_early = exp(dlt),
             p = (1 + sum(abs(pb) >= abs(dlt))) / (1 + NP), n_early = sum(ne), n_late = sum(nl),
             species_early = sum(fe$keep), species_late = sum(fl$keep))
print(lapply(out[c("no_natural", "its", "era")], round, 4)); print(round(apply(out$boot, 2, range), 3))
# named vectors lose their names in write_json: store them as named lists
for (nm in c("main", "no_natural", "its", "era")) out[[nm]] <- as.list(out[[nm]])
write_json(out, sprintf("../data/species_pairs_real%s.json", sfx), digits = 6, auto_unbox = TRUE)
