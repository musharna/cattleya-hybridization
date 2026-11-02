# Figures for pre-registration Y (amendments 1-2): decade path in appearance space, and the selection (S) vs
# transmission (T) split per decade on the mini-catt and trajectory axes, both engines.
# Inputs: data/Y_split.csv, data/Y_decade_pcs.csv (analysis/Y_tables.py). Run from analysis/.
source("house_theme.R")
sp <- read.csv("../data/Y_split.csv")
pc <- read.csv("../data/Y_decade_pcs.csv")
eng <- c(v8 = "orchid-clip-v8", dino = "DINOv2")
sp$engine <- factor(eng[sp$engine], levels = eng)
sp$quantity <- factor(ifelse(sp$quantity == "S", "Selection S (parents' midpoint - pool)",
                             "Transmission T (offspring - parents' midpoint)"))
sp$small <- sp$n < 20

for (e in names(eng)) {
  d <- pc[pc$engine == e, ]
  p <- ggplot(d, aes(PC1, PC2)) +
    geom_path(colour = pal_other) +
    geom_point(aes(size = n), colour = pal_cat[1]) +
    ggrepel::geom_text_repel(aes(label = paste0(decade, "s")), size = 3, seed = 1, min.segment.length = 0.3) +
    scale_size_area(max_size = 5, name = "Grex photos") +
    labs(title = paste0("Decade means of registered Cattleya photos (", eng[[e]], ")"),
         subtitle = "First two principal components of the decade means; one point per registration decade",
         x = "PC1", y = "PC2") +
    theme_house()
  save_house(p, sprintf("../figures/deeper/Y_decade_path_%s.png", e), w = 8, h = 6)
}

split_plot <- function(ax, title, sub) {
  d <- sp[sp$axis == ax, ]
  ggplot(d, aes(decade, estimate, colour = engine)) +
    geom_hline(yintercept = 0, colour = pal_null) +
    geom_pointrange(aes(ymin = lo, ymax = hi, alpha = small), position = position_dodge(width = 4), size = 0.3,
                    show.legend = c(alpha = FALSE)) +
    scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +  # n < 20: faded, unlabelled
    scale_colour_manual(values = unname(pal_cat[1:2]), name = NULL) +
    facet_wrap(~quantity, ncol = 1) +
    labs(title = title, subtitle = sub, x = "Registration decade", y = "Projection on axis (95% parent-unit CI)") +
    theme_house() + theme(legend.position = "bottom", plot.title.position = "plot")
}
save_house(split_plot("minicatt", "Mini-catt axis: the move came from which parents breeders chose",
                      "Axis = held-out photos of 4 founder species minus the labiata group; + = toward the founders"),
           "../figures/deeper/Y_split_minicatt.png", w = 9, h = 7)
save_house(split_plot("trajectory", "Trajectory axis: the two engines disagree on the sign",
                      "Axis = PC1 of held-out decade means; its sign rule (2020s above 1880s) is nearly a tie for v8"),
           "../figures/deeper/Y_split_trajectory.png", w = 9, h = 7)
