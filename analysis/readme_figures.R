# Simple versions of the W and Y figures for the README, drawn from the derived tables only.
# Usage, from analysis/: Rscript readme_figures.R [OUT_DIR]   (default ../figures_regenerated/readme)
# Inputs: ../data/species_pairs_readout.json (W deciles), ../data/Y_split.csv (Y selection, v8, mini-catt axis)
suppressPackageStartupMessages({ library(jsonlite) })
source("theme_cattleya.R")
out <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "../figures_regenerated/readme"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# W: how much more or less often each tenth of species pairs was crossed than expected
w <- fromJSON("../data/species_pairs_readout.json")$figure
stopifnot(nrow(w) == 10, identical(as.integer(w$decile), 1:10))
w$pct <- 100 * (w$obs - 1)
w$dir <- ifelse(w$pct > 0, "more", "fewer")
pw <- ggplot(w, aes(decile, pct, fill = dir)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, colour = ink[["secondary"]]) +
  scale_fill_manual(values = c(more = seq_blue[["dark"]], fewer = ink[["muted"]]), guide = "none") +
  scale_x_continuous(breaks = c(1, 10), labels = c("most closely\nrelated", "least\nrelated")) +
  scale_y_continuous(labels = function(x) ifelse(round(x) == 0, "0%", sprintf("%+d%%", round(x)))) +
  labs(x = "Species pairs, in ten equal groups by relatedness",
       y = "Crosses compared\nwith expected",
       title = "Breeders crossed close relatives more often") +
  theme_cattleya()
save_fig(pw, file.path(out, "W_simple.png"), width = 7, height = 4)

# Y: which way the parents breeders chose leaned, by decade (selection S, orchid-clip-v8, mini-catt axis)
y <- read.csv("../data/Y_split.csv")
y <- y[y$engine == "v8" & y$axis == "minicatt" & y$quantity == "S", ]
stopifnot(nrow(y) == 15, all(diff(y$decade) == 10))
py <- ggplot(y, aes(decade, estimate)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = seq_blue[["light"]], alpha = 0.6) +
  geom_hline(yintercept = 0, colour = ink[["secondary"]]) +
  geom_line(colour = seq_blue[["dark"]], linewidth = 1) +
  geom_point(colour = seq_blue[["dark"]], size = 2.2) +
  scale_y_continuous(breaks = c(-0.08, 0, 0.08),
                     labels = c("labiata\ngroup", "same as all\navailable", "mini-catt\nfounders")) +
  scale_x_continuous(breaks = seq(1880, 2020, 20)) +
  labs(x = "Registration decade", y = "Parents breeders\nchose looked\nmore like",
       title = "Parent choice moved from labiata types to mini-catts") +
  theme_cattleya() + theme(panel.grid.major.y = element_blank())
save_fig(py, file.path(out, "Y_simple.png"), width = 7, height = 4)
cat("wrote", file.path(out, c("W_simple.png", "Y_simple.png")), "\n")
