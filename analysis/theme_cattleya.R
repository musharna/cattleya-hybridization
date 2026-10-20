# House style for every figure in this repository. Plotting scripts source this
# file; restyle here, never inline.
#
# Categorical colours are the first three slots of a palette validated for
# colour-vision deficiency on all pairs (worst CVD ΔE 9.2). The aqua slot sits
# below 3:1 contrast on white, so every series is also labelled directly.

suppressPackageStartupMessages({ library(ggplot2) })

ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984",
         grid = "#e8e7e3", band = "#f1efe9")

hybrid_type_cols <- c(
  "Primary"        = "#2a78d6",   # species × species
  "Complex"        = "#eb6834",   # hybrid × species
  "Highly complex" = "#1baf7a"    # hybrid × hybrid
)

# One-hue ramp for magnitudes (light -> dark blue).
seq_blue <- c(light = "#c9dcf3", mid = "#7fa9e2", dark = "#2a78d6")

corsage_era <- c(start = 1920, end = 1950)

theme_cattleya <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text               = element_text(colour = ink[["primary"]]),
      plot.title         = element_text(face = "bold", size = rel(1.15),
                                        margin = margin(b = 4)),
      plot.title.position = "plot",
      plot.subtitle      = element_text(colour = ink[["secondary"]],
                                        margin = margin(b = 10)),
      plot.caption       = element_text(colour = ink[["muted"]], size = rel(0.75),
                                        hjust = 0),
      plot.caption.position = "plot",
      axis.title         = element_text(colour = ink[["secondary"]], size = rel(0.9)),
      axis.text          = element_text(colour = ink[["secondary"]]),
      panel.grid.major   = element_line(colour = ink[["grid"]], linewidth = 0.35),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x       = element_line(colour = ink[["grid"]], linewidth = 0.35),
      strip.text         = element_text(colour = ink[["primary"]], face = "bold",
                                        hjust = 0),
      legend.position    = "none",
      plot.background    = element_rect(fill = "white", colour = NA),
      plot.margin        = margin(12, 16, 10, 12)
    )
}

# Shaded band + label for the corsage era, for a panel with a year x-axis.
corsage_layers <- function(label_y, label = "corsage era") {
  list(
    annotate("rect", xmin = corsage_era[["start"]], xmax = corsage_era[["end"]],
             ymin = -Inf, ymax = Inf, fill = ink[["band"]]),
    annotate("text", x = mean(corsage_era), y = label_y, label = label,
             colour = ink[["muted"]], size = 3.2, vjust = 1)
  )
}

year_axis <- function() {
  scale_x_continuous(breaks = seq(1850, 2020, by = 25), expand = expansion(mult = c(0.01, 0.02)))
}

save_fig <- function(plot, path, width = 8, height = 4.8) {
  ggsave(path, plot, width = width, height = height, dpi = 200,
         device = ragg::agg_png, bg = "white")
}
