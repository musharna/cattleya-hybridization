# House style for every figure in this repo: source("house_theme.R") from analysis/, then
# add theme_house() and use the palettes below. Restyle here, never inline.
suppressPackageStartupMessages(library(ggplot2))

theme_house <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "right")
}

# categorical palette (Okabe-Ito, colour-blind safe), then neutrals for "other"/"unknown"
pal_cat <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000")
pal_other <- "grey75"
pal_unknown <- "grey45"
pal_obs <- "#D55E00"
pal_null <- "grey60"

save_house <- function(p, path, w = 8, h = 5) ggsave(path, p, width = w, height = h, dpi = 200, bg = "white")
