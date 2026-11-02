#!/usr/bin/env Rscript
# Real-data run for prereg X: main estimate, sensitivities (a) Registrant, (b) walkeriana alone, (c) any adopted
# contact (0/1), and the descriptions (Bass fit, first-10 share, per-species first use), plus one figure.
# Loads breeder_diffusion.R unchanged as a library. Writes aggregates only (no breeder names).
# Usage, from analysis/: Rscript breeder_diffusion_real.R
# Output: ../data/breeder_diffusion_real.json, ../figures_regenerated/deeper/X_breeder_diffusion.png
commandArgs <- function(...) "lib"
source("breeder_diffusion.R")
rm(commandArgs)
source("house_theme.R")

out <- list(main = as.list(fit_soc(d_real)))
sr <- load_set("registrant"); rownames(sr$rows) <- seq_len(nrow(sr$rows))
out$registrant <- as.list(fit_soc(design(sr, sr$A)))
sw <- load_set("originator", "adopt_year_walkeriana"); rownames(sw$rows) <- seq_len(nrow(sw$rows))
out$walkeriana <- as.list(fit_soc(design(sw, sw$A)))
out$binary <- as.list(fit_soc(design(s, s$A, binary = TRUE)))
print(do.call(rbind, lapply(out, function(x) round(unlist(x), 4))))

# descriptions
ad <- read.csv(file.path(BD, "originator", "adoption.csv"), stringsAsFactors = FALSE)
yrs <- 1950:2024
new <- sapply(yrs, function(y) sum(ad$adopt_year == y, na.rm = TRUE))
N0 <- sum(ad$adopt_year < 1950, na.rm = TRUE)
bass <- tryCatch({
  cum <- N0 + c(0, cumsum(new)[-length(new)])
  f <- nls(new ~ (p + q * (cum - N0) / m) * (m - (cum - N0)), start = list(p = 0.01, q = 0.1, m = 500),
           lower = c(1e-6, 0, sum(new)), upper = c(1, 5, 1e5), algorithm = "port")
  as.list(coef(f))
}, error = function(e) list(error = conditionMessage(e)))
post <- ad[!is.na(ad$adopt_year) & ad$adopt_year >= 1950, ]
post <- post[order(post$adopt_year, post$breeder), ]
first10 <- head(post, 10)
out$describe <- list(
  bass = bass, adopters_from1950 = nrow(post),
  first10_share_of_minicatt_grexes = sum(first10$n_mini_from1950) / sum(ad$n_mini_from1950),
  first10_years = range(first10$adopt_year), top10_by_volume_share = sum(head(sort(ad$n_mini_from1950, TRUE), 10)) / sum(ad$n_mini_from1950),
  first_use_by_decade = lapply(c(aclandiae = "adopt_aclandiae", coccinea = "adopt_coccinea", intermedia = "adopt_intermedia",
                                 walkeriana = "adopt_walkeriana"),
                               function(cl) as.list(table(factor(floor(ad[[cl]] / 10) * 10, levels = seq(1850, 2020, 10))))))
str(out$describe[1:5])

# figure: observed adoption rate per breeder-year by adopted contacts, against the asocial model's expectation
ac <- adopted_contacts(s, s$A)[as.integer(rownames(d_real))]
d_real$bin <- cut(ac, c(-Inf, 0, 1, 3, Inf), labels = c("0", "1", "2-3", "4+"))
d_real$expect <- predict(asoc, type = "response")
tab <- do.call(rbind, lapply(split(d_real, d_real$bin), function(g) {
  ci <- prop.test(sum(g$adopt), nrow(g))$conf.int
  data.frame(n = nrow(g), obs = mean(g$adopt), lo = ci[1], hi = ci[2], expect = mean(g$expect))
}))
tab$bin <- factor(rownames(tab), levels = c("0", "1", "2-3", "4+"))
out$figure <- tab
m <- out$main
p <- ggplot(tab, aes(bin)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = pal_cat[1]) +
  geom_point(aes(y = obs), colour = pal_cat[1], size = 3) +
  geom_point(aes(y = expect), shape = 4, size = 3, colour = pal_null) +
  geom_text(aes(y = 0, label = paste0("n = ", n)), vjust = 1.5, size = 3, colour = "grey40") +
  labs(x = "Linked breeders who had already used a mini-catt founder species", y = "Share adopting that year",
       title = "Did breeders take up the mini-catt species after their contacts did?",
       subtitle = sprintf("Dots: observed (95%% CI). Crosses: expected from year, activity and number of contacts alone.\nOdds ratio per log(1 + adopted contacts) %.2f (%.2f-%.2f), p %.2f; %d adoptions 1950-2024",
                          m$or, m$lo, m$hi, m$p, m$n_events)) +
  theme_house() + theme(plot.title.position = "plot")
dir.create("../figures_regenerated/deeper", recursive = TRUE, showWarnings = FALSE)
save_house(p, "../figures_regenerated/deeper/X_breeder_diffusion.png", w = 8, h = 4.5)
write_json(out, "../data/breeder_diffusion_real.json", digits = 6, auto_unbox = TRUE)
