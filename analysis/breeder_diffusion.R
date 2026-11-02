#!/usr/bin/env Rscript
# Did breeders start using the mini-catt founder species sooner when breeders they were linked to already had?
# (prereg/X_breeder_diffusion.md)
#
# What: discrete-time logistic hazard over breeder-years 1950-2024 at risk of first direct use of a founder
# species; exposure = log1p(adopted contacts) with year spline, activity and total contacts as controls.
# Planted truth simulates adoption year by year on the real network and risk set.
# Inputs (local, name people): ~/orchid-data/cattleya_deeper/breeders/{originator,registrant}/
# Usage, from analysis/:
#   Rscript breeder_diffusion.R selftest
#   Rscript breeder_diffusion.R sim ARM START N     # ARM: or1 | or1.5 | or2 | conf0.5 | conf1
#   Rscript breeder_diffusion.R real
# Output: ../data/breeder_diffusion_sim_<ARM>_<START>.json, ../data/breeder_diffusion_real.json (aggregates only)
suppressPackageStartupMessages({library(splines); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "selftest"
BD <- path.expand("~/orchid-data/cattleya_deeper/breeders")

load_set <- function(field = "originator", adopt_col = "adopt_year") {
  by <- read.csv(file.path(BD, field, "breeder_years.csv"), stringsAsFactors = FALSE)
  ed <- read.csv(file.path(BD, field, "edges.csv"), stringsAsFactors = FALSE)
  ad <- read.csv(file.path(BD, field, "adoption.csv"), stringsAsFactors = FALSE)
  br <- ad$breeder; A <- ad[[adopt_col]]; A[is.na(A)] <- Inf
  bi <- match(by$breeder, br); stopifnot(!anyNA(bi))
  # activity: grexes that year and before
  by <- by[order(bi, by$year), ]; bi <- match(by$breeder, br)
  by$n_before <- ave(by$n, bi, FUN = function(x) cumsum(x) - x)
  rows <- by[by$year >= 1950 & by$year <= 2024 & A[bi] >= 1950, ]   # pre-1950 adopters never at risk
  rows$b <- match(rows$breeder, br)
  # contacts of each row's breeder made before that year (long table: row, contact)
  ed$b <- match(ed$breeder, br); ed$c <- match(ed$contact, br); stopifnot(!anyNA(ed$b), !anyNA(ed$c))
  cand <- merge(data.frame(r = seq_len(nrow(rows)), b = rows$b, year = rows$year), ed[, c("b", "c", "first_year")], by = "b")
  long <- cand[cand$first_year < cand$year, c("r", "c", "year")]
  rows$contacts <- tabulate(long$r, nrow(rows))
  # static all-time contact lists, for the confound trait
  list(rows = rows, long = long, A = A, br = br, ed = ed, n_br = length(br))
}

adopted_contacts <- function(s, A) tabulate(s$long$r[A[s$long$c] < s$long$year], nrow(s$rows))

design <- function(s, A, binary = FALSE) {
  d <- s$rows; d <- d[A[d$b] >= d$year, ]                          # still at risk
  ac <- adopted_contacts(s, A)[as.integer(rownames(d))]
  d$adopt <- as.integer(A[d$b] == d$year)
  d$x <- if (binary) as.numeric(ac > 0) else log1p(ac)
  d
}
F_SOC <- adopt ~ ns(year, 4) + log1p(n) + log1p(n_before) + log1p(contacts) + x
F_ASOC <- adopt ~ ns(year, 4) + log1p(n) + log1p(n_before) + log1p(contacts)
fit_soc <- function(d) {
  f <- glm(F_SOC, family = binomial(), data = d)
  co <- summary(f)$coefficients["x", ]
  c(beta = co[[1]], se = co[[2]], or = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]),
    hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]], n_events = sum(d$adopt), n_rows = nrow(d))
}

s <- load_set()
stopifnot(all(as.integer(rownames(s$rows)) >= 1))
rownames(s$rows) <- seq_len(nrow(s$rows))
d_real <- design(s, s$A)
asoc <- glm(F_ASOC, family = binomial(), data = d_real)          # baseline for the simulations only

# simulate adoption year by year; pre-1950 adopters keep their real year
sim_A <- function(s, log_or, gamma = 0, z = NULL) {
  A <- ifelse(s$A < 1950, s$A, Inf)
  r <- s$rows
  lp0 <- predict(asoc, newdata = r, type = "link")
  for (y in 1950:2024) {
    at <- which(r$year == y & A[r$b] > y - 0.5)
    if (!length(at)) next
    ac <- tabulate(s$long$r[s$long$r %in% at & A[s$long$c] < y], nrow(r))[at]
    lp <- lp0[at] + log_or * log1p(ac) + if (gamma != 0) gamma * z[r$b[at]] else 0
    hit <- at[runif(length(at)) < plogis(lp)]
    A[r$b[hit]] <- y
  }
  A
}
conf_trait <- function(s) {
  e <- rnorm(s$n_br)
  m <- tapply(e[s$ed$c], s$ed$b, mean); nb <- numeric(s$n_br); nb[as.integer(names(m))] <- m
  z <- e + nb; (z - mean(z)) / sd(z)
}

if (mode == "selftest") {
  cat(sprintf("risk rows %d, breeders at risk %d, events %d, long contacts %d\n", nrow(d_real), length(unique(d_real$b)),
              sum(d_real$adopt), nrow(s$long)))
  stopifnot(sum(d_real$adopt) == sum(is.finite(s$A) & s$A >= 1950 & s$A <= 2024 & s$A %in% s$rows$year))
  set.seed(1); t0 <- proc.time()[["elapsed"]]
  A2 <- sim_A(s, log(2)); r2 <- fit_soc(design(s, A2))
  A1 <- sim_A(s, 0); r1 <- fit_soc(design(s, A1))
  cat(sprintf("planted OR 2: recovered %.2f (%.2f-%.2f), events %d; planted OR 1: %.2f p %.3f; %.1f s for 2 sims\n",
              r2[["or"]], r2[["lo"]], r2[["hi"]], r2[["n_events"]], r1[["or"]], r1[["p"]], proc.time()[["elapsed"]] - t0))
}
if (mode == "sim") {
  arm <- args[2]; st <- as.integer(args[3]); nn <- as.integer(args[4])
  spec <- list(or1 = c(0, 0, 6000), or1.5 = c(log(1.5), 0, 7000), or2 = c(log(2), 0, 8000),
               conf0.5 = c(0, 0.5, 9000), conf1 = c(0, 1, 9500))[[arm]]
  stopifnot(!is.null(spec))
  res <- t(sapply(st:(st + nn - 1), function(k) {
    set.seed(spec[3] + k)
    z <- if (spec[2] != 0) conf_trait(s) else NULL
    out <- fit_soc(design(s, sim_A(s, spec[1], spec[2], z))); cat(k, round(out, 4), "\n"); out
  }))
  write_json(list(arm = arm, start = st, reps = res), sprintf("../data/breeder_diffusion_sim_%s_%d.json", arm, st), digits = 6)
}
