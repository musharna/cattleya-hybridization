# Adapter for the later-analysis scripts: they call theme_house(), pal_* and save_house(), which now
# take their look from the repository's one house style in theme_cattleya.R. Restyle there, not here.
source("theme_cattleya.R")

# Same as theme_cattleya() but keeps legends, which several of these figures need.
theme_house <- function(base_size = 12) theme_cattleya(base_size) + theme(legend.position = "right")

# Categorical: the three house colours first, then four more Okabe-Ito hues for figures with up to 8 series.
pal_cat <- c(unname(hybrid_type_cols), "#CC79A7", "#56B4E9", "#F0E442", "#000000", "#8C6D31")
stopifnot(length(pal_cat) == 8, !anyDuplicated(pal_cat))
pal_other <- "#c9c8c3"
pal_unknown <- ink[["secondary"]]
pal_obs <- hybrid_type_cols[["Complex"]]
pal_null <- ink[["muted"]]

save_house <- function(p, path, w = 8, h = 5) save_fig(p, path, width = w, height = h)
