#!/usr/bin/env Rscript
# Regenerate the hybridization-trend figures from the derived tables in data/.
# Needs no access to the raw orchidroots scrape.
#
#   Rscript figures_from_derived.R [output-dir]
#
# Covers Figure 1, Figure 2, Figure S1 and Figure S2. Figure 3 and Table 1
# come from species_tree_signal.R.

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

out <- if (length(commandArgs(trailingOnly = TRUE)) >= 1)
  commandArgs(trailingOnly = TRUE)[1] else "../figures_regenerated"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

trend    <- read.csv("../data/hybrids_by_year_type.csv")
ancestors<- read.csv("../data/ancestors_by_year.csv")
per_year <- read.csv("../data/hybrids_per_year.csv")

# --- Figure 1: hybridization trends by type, with the Corsage Era band ------
trend$HybridTypeLabel <- dplyr::recode(trend$HybridType,
  "Primary Hybrid"        = "Primary (species × species)",
  "Complex Hybrid"        = "Complex (hybrid × species)",
  "Highly Complex Hybrid" = "Highly Complex (hybrid × hybrid)")

corsage_era <- data.frame(xmin = 1920, xmax = 1950, ymin = -Inf, ymax = Inf,
                          era = "Cattleya Corsage Era")

p1 <- ggplot() +
  geom_rect(data = corsage_era,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = era),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = trend,
            aes(x = Year, y = Count, color = HybridTypeLabel, group = HybridTypeLabel),
            linewidth = 1) +
  geom_point(data = trend, aes(x = Year, y = Count, color = HybridTypeLabel)) +
  scale_fill_manual(name = "", values = c("Cattleya Corsage Era" = "purple")) +
  labs(title = "Hybridization Trends Over Time", x = "Year", y = "Count",
       color = "Hybrid Type") +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(alpha = 0.2)))
ggsave(file.path(out, "fig1_hybridization_trends_by_type.png"), p1,
       width = 9, height = 5.5, dpi = 150)

# --- Figure 2: ancestor count distribution by year -------------------------
p2 <- ggplot(ancestors, aes(x = as.factor(Year), y = Ancestors)) +
  geom_boxplot(fill = "lightblue", outlier.size = 0.8) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "red", linewidth = 1, linetype = "dashed") +
  labs(title = "Number of Ancestors by Year", x = "Year", y = "Number of Ancestors") +
  scale_x_discrete(breaks = seq(1850, max(ancestors$Year), by = 10)) +
  theme_minimal()
ggsave(file.path(out, "fig2_ancestors_by_year.png"), p2,
       width = 9, height = 5.5, dpi = 150)

# --- Figure S1: hybrids registered per year --------------------------------
pS1 <- ggplot(per_year, aes(x = Year, y = Hybrids)) +
  geom_line(color = "blue") + geom_point(color = "red") +
  labs(title = "Number of Cattleya Hybrids per Year", x = "Year", y = "Number of Hybrids") +
  scale_x_continuous(breaks = seq(1850, max(per_year$Year), by = 10)) +
  theme_minimal()
ggsave(file.path(out, "figS1_hybrids_per_year.png"), pS1,
       width = 9, height = 5.5, dpi = 150)

# --- Figure S2: cumulative registered hybrids ------------------------------
pS2 <- ggplot(per_year, aes(x = Year, y = Cumulative)) +
  geom_line(color = "steelblue") +
  labs(title = "Cumulative Number of Registered Hybrids Over Time",
       x = "Year", y = "Cumulative Hybrid Count") +
  theme_minimal()
ggsave(file.path(out, "figS2_cumulative_hybrids.png"), pS2,
       width = 11, height = 6, dpi = 150)

cat("wrote 4 figures to", normalizePath(out), "\n")
