# Prereg C Amendment 8: power curve for test 3 on the combined molecular tree.
# Logistic regression of detect (H2 95% lower bound >= 0.05) on planted lambda; lambda80 = lambda at fitted p 0.8,
# with a 2,000-resample bootstrap within each lambda. Unconverged fits (Rhat >= 1.01) are listed and left out.
#   Rscript power_curve.R            # real: e5 sims, e7 shards replacing theirs, e8 shards added
#   Rscript power_curve.R selftest   # planted logistic with a known lambda80; must be recovered
# Output: ../data/signal_corrected_power_curve.json
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(TRUE)

lambda80 <- function(df) {
  f <- suppressWarnings(glm(detect ~ planted, family = binomial, data = df))
  b <- coef(f)
  if (!is.finite(b[2]) || b[2] <= 0) return(NA_real_)
  (qlogis(0.8) - b[1]) / b[2]
}
boot80 <- function(df, n = 2000) {
  set.seed(1)
  idx <- split(seq_len(nrow(df)), df$planted)
  v <- replicate(n, lambda80(df[unlist(lapply(idx, function(i) i[sample.int(length(i), replace = TRUE)])), ]))
  bad <- is.na(v) | v < 0 | v > 1
  v[bad] <- 1
  list(lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)), n_bad = sum(bad))
}

if (length(args) && args[1] == "selftest") {
  # Known answer: p = plogis(a + b*lambda) with lambda80 = 0.5; 200 fits per lambda so the estimate must land close
  set.seed(7); a <- -3; b <- (qlogis(0.8) - a) / 0.5
  grid <- c(0, 0.3, 0.45, 0.6, 0.8)
  df <- data.frame(planted = rep(grid, each = 200))
  df$detect <- rbinom(nrow(df), 1, plogis(a + b * df$planted))
  est <- lambda80(df); ci <- boot80(df, 500)
  cat(sprintf("selftest: true 0.500, est %.3f [%.3f, %.3f]\n", est, ci$lo, ci$hi))
  stopifnot(abs(est - 0.5) < 0.03, ci$lo < 0.5, ci$hi > 0.5)
  # A flat curve (no dependence on lambda) must not yield an in-range lambda80
  df$detect <- rbinom(nrow(df), 1, 0.3)
  flat <- lambda80(df)
  cat(sprintf("selftest flat: est %s\n", format(flat)))
  stopifnot(is.na(flat) || flat > 1 || flat < 0)
  cat("SELFTEST PASS\n"); quit(status = 0)
}

key <- function(f) sub(".*_(l[0-9.]+_s[0-9]+)\\.json$", "\\1", f)
sf <- Sys.glob("../data/signal_corrected_sim3_mol_e5_l*_s*.json")
s7 <- Sys.glob("../data/signal_corrected_sim3_mol_e7_l*_s*.json")
s8 <- Sys.glob("../data/signal_corrected_sim3_mol_e8_l*_s*.json")
stopifnot(length(sf) == 6, length(s7) == 3, length(s8) == 4, key(s7) %in% key(sf), !key(s8) %in% key(sf))
sf[match(key(s7), key(sf))] <- s7
df <- do.call(rbind, lapply(c(sf, s8), function(f) transform(fromJSON(f)$reps, file = basename(f))))
stopifnot(nrow(df) == 50, all(table(df$planted) == 10))
unconv <- df[df$rhat_max >= 1.01, ]
df <- df[df$rhat_max < 1.01, ]
df$detect <- as.integer(df$lo >= 0.05)
per <- do.call(rbind, lapply(split(df, df$planted), function(d)
  data.frame(planted = d$planted[1], detect = sum(d$detect), n = nrow(d), median_h2 = median(d$median))))
est <- lambda80(df); ci <- boot80(df)
print(per); cat(sprintf("lambda80 %.3f [%.3f, %.3f], bootstrap resamples set to 1: %d; unconverged left out: %d\n",
                        est, ci$lo, ci$hi, ci$n_bad, nrow(unconv)))
write_json(list(per_lambda = per, lambda80 = est, lo = ci$lo, hi = ci$hi, n_boot_bad = ci$n_bad,
                unconverged = unconv), "../data/signal_corrected_power_curve.json", digits = 6, auto_unbox = TRUE)

# Figure: detection rate per planted lambda with the fitted curve and lambda80
source("house_theme.R")
fit <- glm(detect ~ planted, family = binomial, data = df)
grid <- data.frame(planted = seq(0, 0.8, by = 0.01)); grid$p <- predict(fit, grid, type = "response")
per$rate <- per$detect / per$n
p <- ggplot() +
  annotate("rect", xmin = ci$lo, xmax = ci$hi, ymin = -Inf, ymax = Inf, fill = pal_other, alpha = 0.4) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = pal_null) +
  geom_line(data = grid, aes(planted, p), colour = pal_cat[1]) +
  geom_point(data = per, aes(planted, rate), size = 3, colour = pal_cat[1]) +
  geom_text(data = per, aes(planted, rate, label = sprintf("%d/%d", detect, n)), vjust = -1, size = 3.2) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  labs(x = "Planted phylogenetic signal (lambda)", y = "Share of simulations detecting signal",
       title = "How strong clustering must be for the test to see it",
       subtitle = sprintf("Detected 4 times in 5 at lambda %.2f (95%% interval %.2f-%.2f, shaded)", est, ci$lo, ci$hi)) +
  theme_house()
dir.create("../figures_regenerated/deeper", recursive = TRUE, showWarnings = FALSE)
save_house(p, "../figures_regenerated/deeper/C_power_curve.png", w = 7, h = 4.5)
