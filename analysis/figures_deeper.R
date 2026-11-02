#!/usr/bin/env Rscript
# Figures for the "deeper" analyses (G, A, E, B/D simulations, C), one graph per file.
# Usage, from analysis/: Rscript figures_deeper.R G A E sims C   (any subset; a missing input fails)
# Output: ../figures_regenerated/deeper/
suppressPackageStartupMessages({library(jsonlite); library(dplyr); library(tidyr)})
source("house_theme.R")
args <- commandArgs(TRUE)
out <- "../figures_regenerated/deeper"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if ("G" %in% args) {
  g <- read.csv("../data/qc_summary.csv")
  p <- ggplot(g, aes(reorder(class, records_flagged), records_flagged)) +
    geom_col(fill = pal_cat[1]) + coord_flip() +
    labs(x = NULL, y = "Records flagged (of 19,577)", title = "Register QC flags by class") + theme_house()
  save_house(p, file.path(out, "G_qc_flags.png"), h = 4)
}

if ("A" %in% args) {
  sh <- read.csv("../data/founder_shares_by_decade.csv")
  fd <- read.csv("../data/founders_by_decade.csv")
  keep <- fd$decade[fd$n >= 20]
  sh <- filter(sh, decade %in% keep)
  top <- sh %>% filter(species != "unknown") %>% group_by(species) %>% summarise(s = sum(share)) %>%
    arrange(desc(s)) %>% head(8) %>% pull(species)
  sh2 <- sh %>% mutate(group = ifelse(species == "unknown", "unknown", ifelse(species %in% top, sub("Cattleya ", "C. ", species), "other species"))) %>%
    group_by(decade, group) %>% summarise(share = sum(share), .groups = "drop") %>%
    complete(decade, group, fill = list(share = 0))  # an absent species is share 0, not a gap
  lv <- c(sub("Cattleya ", "C. ", top), "other species", "unknown")
  sh2$group <- factor(sh2$group, levels = lv)
  p <- ggplot(sh2, aes(decade, share, fill = group)) + geom_area() +
    scale_fill_manual(values = setNames(c(pal_cat, pal_other, pal_unknown), lv), name = NULL) +
    labs(x = "Registration decade", y = "Expected founder-species share",
         title = "Founder-species composition of each decade's registrations") + theme_house()
  save_house(p, file.path(out, "A_founder_shares.png"))

  gl <- fd %>% filter(n >= 20) %>% select(decade, mini, labiata) %>%
    pivot_longer(-decade, names_to = "group", values_to = "share") %>%
    mutate(group = recode(group, mini = "mini-catt founders (coccinea, walkeriana, aclandiae, intermedia)",
                          labiata = "labiata group (labiata, mossiae, trianae, warscewiczii, dowiana)"))
  p <- ggplot(gl, aes(decade, share, colour = group)) + geom_line(linewidth = 1) + geom_point() +
    scale_colour_manual(values = pal_cat[c(2, 1)], name = NULL) +
    labs(x = "Registration decade", y = "Cohort share", title = "Positive control: mini-catt rise, labiata-group fall") +
    theme_house() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_house(p, file.path(out, "A_control_minicatt.png"))

  res <- fromJSON("../data/A_result.json")
  band <- bind_rows(lapply(names(res$null$band), function(d) {
    b <- res$null$band[[d]]
    bind_rows(lapply(names(b), function(m) data.frame(decade = as.integer(d), metric = m, obs = b[[m]]$obs, lo = b[[m]]$lo, hi = b[[m]]$hi)))
  }))
  lab <- c(f_e = "Effective number of founder species f_e", ratio = "f_e / f_a (frozen bottleneck test)",
           top_grex_contrib = "Largest single-grex contribution")
  ttl <- c(f_e = "Founder diversity f_e", ratio = "Bottleneck ratio f_e / f_a", top_grex_contrib = "Top grex contribution")
  for (m in names(lab)) {
    d <- filter(band, metric == m)
    p <- ggplot(d, aes(decade)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = pal_null, alpha = 0.5) +
      geom_line(aes(y = obs), colour = pal_obs, linewidth = 1) + geom_point(aes(y = obs), colour = pal_obs) +
      labs(x = "Registration decade", y = lab[[m]],
           title = paste0(ttl[[m]], ": register vs edge-swap null"),
           subtitle = sprintf("Line = register; band = 95%% of %d degree-preserving swap draws", res$null$draws)) +
      theme_house()
    save_house(p, file.path(out, sprintf("A_null_%s.png", m)))
  }
}

if ("E" %in% args) {
  e <- fromJSON("../data/e_gate.json")
  t <- bind_rows(lapply(names(e$fisher), function(s) {
    m <- e$fisher[[s]]$table  # fromJSON gives the row-major matrix: rows notho yes/no, cols primary yes/no
    stopifnot(is.matrix(m), all(dim(m) == 2))
    data.frame(stratum = s, primary = c("crossed (registered primary)", "never crossed"),
               share_notho = c(m[1, 1] / (m[1, 1] + m[2, 1]), m[1, 2] / (m[1, 2] + m[2, 2])),
               n = c(m[1, 1] + m[2, 1], m[1, 2] + m[2, 2]), OR = e$fisher[[s]]$odds_ratio, p = e$fisher[[s]]$p)
  }))
  t$stratum <- factor(t$stratum, levels = c("pooled", "same_rank", "different_rank"),
                      labels = c("all pairs", "same subgeneric rank", "different rank"))
  p <- ggplot(t, aes(primary, share_notho, fill = primary)) + geom_col() +
    geom_text(aes(label = sprintf("%.1f%% of %d pairs", 100 * share_notho, n)), vjust = -0.3, size = 3) +
    facet_wrap(~stratum) + scale_fill_manual(values = pal_cat[c(1, 5)], guide = "none") +
    scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
    labs(x = NULL, y = "Share of species pairs with a WCVP nothospecies",
         title = "Pairs that breeders have crossed are more often natural hybrids",
         subtitle = with(t[t$stratum == "all pairs", ][1, ], sprintf("All pairs: odds ratio %.1f, Fisher p = %.1g", OR, p))) +
    theme_house()
  save_house(p, file.path(out, "E_2x2_share.png"), w = 9)
  k <- e$readout1_breakdown
  br <- data.frame(kind = c("independent registered primary", "capitalised remake under the nothospecies epithet",
                            "species-style register entry only", "no register record"),
                   n = c(k$independent, k$remake_under_notho_epithet, k$species_style_only, k$none))
  br$kind <- factor(br$kind, levels = rev(br$kind))
  p <- ggplot(br, aes(kind, n)) + geom_col(fill = pal_cat[1]) + coord_flip() +
    labs(x = NULL, y = sprintf("Nothospecies pairs in the gate universe (n = %d)", sum(br$n)),
         title = "E gate readout 1: register evidence per natural hybrid",
         subtitle = sprintf("Bar >= 75%% independent; observed %.0f%% (%.0f%% with remakes)",
                            100 * br$n[1] / sum(br$n), 100 * sum(br$n[1:2]) / sum(br$n))) + theme_house() + theme(plot.title.position = "plot")
  save_house(p, file.path(out, "E_readout1_breakdown.png"), w = 11, h = 3.5)
}

if ("sims" %in% args) {
  s <- fromJSON("../data/prereg_sims.json")
  b <- s$B$cells %>% mutate(cell = sprintf("a=%d b1=%.2f c=%+.3f", a, b1, c))
  p <- ggplot(b, aes(factor(delta), pass_rate, group = cell, colour = factor(c))) + geom_line() + geom_point() +
    geom_hline(yintercept = c(0.10, 0.80), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = pal_cat[1:3], name = "curvature c") +
    labs(x = "Planted slope change at 1960 (years per year)", y = "Pass rate (200 replicates)",
         title = "B: power and false-positive rate of the 1960 break test",
         subtitle = "Dashed: frozen bars (FPR <= 0.10 at 0, power >= 0.80 at 0.2)") + theme_house()
  save_house(p, file.path(out, "B_sim_power.png"))
  d <- s$D$cells %>% mutate(OR = round(exp(gamma), 2))
  p <- ggplot(d, aes(factor(OR), rate, colour = factor(kappa), group = factor(kappa))) + geom_line() + geom_point() +
    geom_hline(yintercept = c(0.07, 0.80), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = pal_cat[1:2], name = "popularity effect") +
    labs(x = "Planted odds ratio, different vs same rank", y = "Rejection rate (500 replicates)",
         title = "D: power of the fertility-ladder test",
         subtitle = "Dashed: frozen bars (FPR <= 0.07 at OR 1, power >= 0.80 at OR 0.7)") + theme_house()
  save_house(p, file.path(out, "D_sim_power.png"))
}

if ("C" %in% args) {
  # Amendment 5: every C number below uses the corrected exposure (tags *_e5; Amendment 6 refits *_e6). The v1-tree-only figure
  # (C_h2.png) is retired: its planted truth was not rerun, and the v1 tree is a row of the figure below.
  r <- fromJSON("../data/signal_corrected_real_e5.json")
  reg <- function(what, v)  # h2(): median, lo, hi, rhat
    data.frame(what = paste0(what, if (length(v) > 3 && v[4] >= 1.01) " (Rhat > 1.01)" else ""),
               median = v[1], lo = v[2], hi = v[3])

  # Amendment 2: the same test on the molecular (GenBank) trees
  q <- function(x) unname(quantile(x, c(0.5, 0.025, 0.975)))
  rm_ <- function(tag, e = "e5") fromJSON(sprintf("../data/signal_corrected_real_%s_%s.json", tag, e))$test3$with_transferred
  # Amendment 6: the five unconverged bootstrap trees are replaced by their e6 refits
  per_tree <- function(pat) {
    out <- list()
    for (f in Sys.glob(pat)) {
      j <- fromJSON(f, simplifyVector = FALSE)
      for (k in seq_along(j$trees)) out[[as.character(j$trees[[k]])]] <- list(h = unlist(j$h[[k]]), rhat = j$rhat_max[[k]])
    }
    out
  }
  bt <- per_tree("../data/signal_corrected_boot_mol_e5_*.json"); b6 <- per_tree("../data/signal_corrected_boot_mol_e6_*.json")
  stopifnot(length(bt) == 20, length(b6) == 5, names(b6) %in% names(bt))
  bt[names(b6)] <- b6
  boot <- unlist(lapply(bt, `[[`, "h"))
  n_unconv <- sum(sapply(bt, `[[`, "rhat") >= 1.01)
  # Amendment 7: the three rerun shards (*_e7) replace their e5 versions whole
  sf <- Sys.glob("../data/signal_corrected_sim3_mol_e5_l*_s*.json"); s7 <- Sys.glob("../data/signal_corrected_sim3_mol_e7_l*_s*.json")
  key <- function(f) sub(".*_(l[0-9.]+_s[0-9]+)\\.json$", "\\1", f)
  stopifnot(length(sf) == 6, length(s7) == 3, key(s7) %in% key(sf))
  sf[match(key(s7), key(sf))] <- s7
  sm <- bind_rows(lapply(sf, function(f) fromJSON(f)$reps))
  sm <- sm[order(sm$planted), ]
  hm <- rbind(reg("combined tree (91 spp.)", rm_("mol")), reg(sprintf("20 bootstrap trees, pooled%s", if (n_unconv) sprintf(" (%d with Rhat > 1.01)", n_unconv) else ""), q(boot)),
              reg("ITS-only tree", rm_("mol_its", "e6")), reg("plastid-only tree (51 spp.)", rm_("mol_plastid")),
              reg("100 spp., 9 grafted", rm_("mol_all100", "e6")), reg("v1 subgeneric tree", r$test3$with_transferred),
              data.frame(what = sprintf("planted lambda %s%s", sm$planted, ifelse(sm$rhat_max >= 1.01, " (Rhat > 1.01)", "")),
                         median = sm$median, lo = sm$lo, hi = sm$hi))
  hm$id <- seq_len(nrow(hm)); hm$reg <- !grepl("planted", hm$what)
  p <- ggplot(hm, aes(reorder(paste(what, id), -id), median, colour = reg)) +
    geom_pointrange(aes(ymin = lo, ymax = hi)) + coord_flip() +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c(`TRUE` = pal_obs, `FALSE` = pal_null), guide = "none") +
    scale_x_discrete(labels = function(x) sub(" [0-9]+$", "", x)) +
    labs(x = NULL, y = "Phylogenetic H2 (posterior median, 95% CrI)",
         title = "C on a molecular tree: still no supported signal",
         subtitle = "Orange: register, by tree. Grey: planted truth, combined tree. Dashed: bar 0.05") +
    theme_house() + theme(plot.title.position = "plot")
  save_house(p, file.path(out, "C_h2_molecular.png"), h = 7)
}
