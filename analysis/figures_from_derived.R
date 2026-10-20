#!/usr/bin/env Rscript
# Regenerate the hybridization-trend figures from the derived tables in data/.
# Needs no access to the raw orchidroots scrape.
#
#   Rscript figures_from_derived.R [output-dir]
#
# Covers Figure 1, Figure 2, Figure S1 and Figure S2. Figure 3 and Table 1
# come from species_tree_signal.R. Styling lives in theme_cattleya.R.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})
source("theme_cattleya.R")

out <- if (length(commandArgs(trailingOnly = TRUE)) >= 1)
  commandArgs(trailingOnly = TRUE)[1] else "../figures_regenerated"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

trend     <- read.csv("../data/hybrids_by_year_type.csv")
ancestors <- read.csv("../data/ancestors_by_year.csv")
per_year  <- read.csv("../data/hybrids_per_year.csv")

type_labels <- c("Primary Hybrid" = "Primary", "Complex Hybrid" = "Complex",
                 "Highly Complex Hybrid" = "Highly complex")
type_long <- c("Primary" = "Primary\n(species × species)",
               "Complex" = "Complex\n(hybrid × species)",
               "Highly complex" = "Highly complex\n(hybrid × hybrid)")
years <- seq(min(trend$Year), max(trend$Year))

# Years with no registrations are absent from the tables; they are zeros.
trend <- trend |>
  mutate(Type = factor(type_labels[HybridType], levels = names(hybrid_type_cols))) |>
  select(Year, Type, Count) |>
  complete(Year = years, Type, fill = list(Count = 0)) |>
  arrange(Type, Year) |>
  group_by(Type) |>
  mutate(Smooth = as.numeric(stats::filter(Count, rep(1 / 5, 5), sides = 2))) |>
  ungroup()

bin5 <- function(y) 5 * floor(y / 5)

# --- Figure 1: registrations by type, and each type's share ----------------
ends <- trend |> filter(!is.na(Smooth)) |> group_by(Type) |> slice_max(Year, n = 1)

p1a <- ggplot(filter(trend, !is.na(Smooth)), aes(Year, Smooth, colour = Type)) +
  corsage_layers(label_y = 265) +
  geom_line(linewidth = 0.9) +
  geom_text(data = ends, aes(label = type_long[as.character(Type)]),
            hjust = 0, nudge_x = 2, size = 3.1, lineheight = 0.9,
            colour = ink[["primary"]]) +
  scale_colour_manual(values = hybrid_type_cols) +
  year_axis() +
  coord_cartesian(xlim = c(1850, 2048), clip = "off") +
  labs(title = "Hybrid × hybrid crosses took over and peaked in the 1960s",
       subtitle = "Cattleya hybrids registered per year, by kind of cross (5-year moving average)",
       x = NULL, y = "Hybrids per year") +
  theme_cattleya()

share <- trend |>
  mutate(Bin = bin5(Year)) |>
  group_by(Bin, Type) |> summarise(n = sum(Count), .groups = "drop") |>
  group_by(Bin) |> mutate(total = sum(n), Share = n / total) |> ungroup() |>
  filter(total >= 20)

# Label each band at a period where it is wide and away from the panel edge.
label_bin <- c("Primary" = 1890, "Complex" = 1990, "Highly complex" = 1955)
share_labels <- share |>
  arrange(Bin, desc(Type)) |> group_by(Bin) |>
  mutate(y = cumsum(Share) - Share / 2) |> ungroup() |>
  filter(Bin == label_bin[as.character(Type)])

p1b <- ggplot(share, aes(Bin + 2.5, Share, fill = Type)) +
  geom_area(colour = "white", linewidth = 0.5) +
  geom_text(data = share_labels, aes(x = Bin + 2.5, y = y, label = Type), colour = "white",
            fontface = "bold", size = 3.2) +
  scale_fill_manual(values = hybrid_type_cols) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  year_axis() +
  coord_cartesian(xlim = c(1850, 2048)) +
  labs(subtitle = "Share of each 5-year period's registrations (periods with at least 20)",
       x = NULL, y = "Share") +
  theme_cattleya()

save_fig(p1a / p1b + plot_layout(heights = c(3, 2)),
         file.path(out, "fig1_hybridization_trends_by_type.png"), width = 8.5, height = 7)

# --- Figure 2: ancestors per hybrid, 5-year bins ----------------------------
anc <- ancestors |>
  mutate(Bin = bin5(Year)) |>
  group_by(Bin) |>
  summarise(n = n(), mean = mean(Ancestors), median = median(Ancestors),
            q10 = quantile(Ancestors, 0.10), q25 = quantile(Ancestors, 0.25),
            q75 = quantile(Ancestors, 0.75), q90 = quantile(Ancestors, 0.90)) |>
  filter(n >= 20) |>
  mutate(x = Bin + 2.5)
last <- slice_max(anc, x, n = 1)

p2 <- ggplot(anc, aes(x)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = seq_blue[["light"]], alpha = 0.6) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = seq_blue[["mid"]], alpha = 0.6) +
  geom_line(aes(y = median), colour = seq_blue[["dark"]], linewidth = 0.9) +
  geom_line(aes(y = mean), colour = ink[["primary"]], linewidth = 0.6, linetype = "22") +
  annotate("text", x = last$x + 2, y = c(last$q90, last$q75, last$mean, last$median),
           label = c("90th percentile", "75th percentile", "mean", "median"),
           hjust = 0, size = 3.1,
           colour = c(ink[["secondary"]], ink[["secondary"]], ink[["primary"]], ink[["primary"]])) +
  year_axis() +
  coord_cartesian(xlim = c(1850, 2045), clip = "off") +
  labs(title = "Hybrids grew more complex until the 1970s, then levelled off",
       subtitle = "Ancestors per registered Cattleya hybrid, 5-year bins (shaded: middle 50% and middle 80%)",
       x = NULL, y = "Ancestors per hybrid") +
  theme_cattleya()
save_fig(p2, file.path(out, "fig2_ancestors_by_year.png"), width = 8.5, height = 4.8)

# --- Figure S1: hybrids registered per year --------------------------------
py <- tibble(Year = years) |> left_join(per_year, by = "Year") |>
  mutate(Hybrids = coalesce(as.numeric(Hybrids), 0))
peak <- slice_max(py, Hybrids, n = 1)

pS1 <- ggplot(py, aes(Year, Hybrids)) +
  corsage_layers(label_y = 330) +
  geom_line(colour = seq_blue[["dark"]], linewidth = 0.6) +
  annotate("point", x = peak$Year, y = peak$Hybrids, colour = seq_blue[["dark"]], size = 2.5) +
  annotate("text", x = peak$Year + 3, y = peak$Hybrids,
           label = sprintf("%d: %d", peak$Year, peak$Hybrids), hjust = 0, size = 3.1,
           colour = ink[["secondary"]]) +
  year_axis() +
  labs(title = "Registrations peaked around 1960",
       subtitle = "Cattleya hybrids registered per year", x = NULL, y = "Hybrids per year") +
  theme_cattleya()
save_fig(pS1, file.path(out, "figS1_hybrids_per_year.png"))

# --- Figure S2: cumulative registered hybrids ------------------------------
endS2 <- slice_max(per_year, Year, n = 1)
pS2 <- ggplot(per_year, aes(Year, Cumulative)) +
  geom_line(colour = seq_blue[["dark"]], linewidth = 0.9) +
  annotate("text", x = endS2$Year - 2, y = endS2$Cumulative,
           label = format(endS2$Cumulative, big.mark = ","), hjust = 1, vjust = -0.6,
           size = 3.1, colour = ink[["secondary"]]) +
  scale_y_continuous(labels = scales::comma) +
  year_axis() +
  labs(title = "Registered Cattleya hybrids over time",
       subtitle = "Running total of hybrids with a registration year", x = NULL,
       y = "Registered hybrids") +
  theme_cattleya()
save_fig(pS2, file.path(out, "figS2_cumulative_hybrids.png"))

cat("wrote 4 figures to", normalizePath(out), "\n")
