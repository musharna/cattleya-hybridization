#!/usr/bin/env Rscript
# Phylogenetic signal tests on the subgeneric-rank and species-level trees.
#
# This reproduces Figure 3 and both rows of Table 1 from the derived tables in
# data/ alone -- no access to the raw scrape is needed.
#
# The procedure is taken verbatim from the submitted analysis. The only change is the input: that file read the raw species CSV and
# computed `LowestRank` from the multi-line "Subgeneric Ranks" field, whereas
# this script reads the already-reduced `LastSubgenus` column, which
# derive_tables.R computes the same way (last non-empty line of that field).
#
# Usage, from analysis/:   Rscript species_tree_signal.R

suppressPackageStartupMessages({
  library(ape); library(phytools); library(phylosignal)
  library(phylobase); library(ggplot2); library(patchwork)
})
source("theme_cattleya.R")

set.seed(1)                       # permutation tests below use 999 replicates
out_dir <- "../figures_regenerated"
dir.create(out_dir, showWarnings = FALSE)

# ---- inputs ---------------------------------------------------------------
tree_lowestrank <- read.tree("../data/cattleya_subgeneric.nwk")
sp <- read.csv("../data/species_hybrid_counts.csv", stringsAsFactors = FALSE)
sp$SpeciesName <- gsub(" ", "_", sp$Binomial)

cat(sprintf("subgeneric tree: %d tips\n", length(tree_lowestrank$tip.label)))
cat(sprintf("species records: %d\n", nrow(sp)))

# ---- build the species-level tree -----------------------------------------
# Each subgeneric rank tip is replaced by a star polytomy of its species.
rank_to_species <- split(sp$SpeciesName, sp$LastSubgenus)

# A rank holding a single species becomes that species' tip: the submitted
# code expressed this as bind.tree() of a zero-edge phylo, which under
# ape 5.7.1 does not terminate. Renaming the tip is what that bind resolved
# to -- binding the first singleton left the tip count unchanged at 10,
# i.e. it substituted the tip rather than adding one.
tree_species <- tree_lowestrank
for (rank in names(rank_to_species)) {
  s <- rank_to_species[[rank]]
  idx <- which(tree_species$tip.label == rank)
  if (length(idx) == 0) next          # rank absent from the subgeneric tree
  if (length(s) == 1) {
    tree_species$tip.label[idx] <- s
  } else {
    tree_species <- bind.tree(
      tree_species, stree(length(s), tip.label = s, type = "star"),
      where = idx)
  }
}
tree_species <- drop.tip(
  tree_species,
  intersect(tree_species$tip.label, tree_lowestrank$tip.label))

# Ranks absent from the subgeneric tree (Crispae, Stellata) never get bound,
# so their species drop out here -- matching the submitted analysis, which
# excluded Stellata explicitly and restricted to tips present in the tree.
dropped <- setdiff(sp$SpeciesName, tree_species$tip.label)
cat(sprintf("species-level tree: %d tips (%d dropped: ranks not in tree)\n",
            length(tree_species$tip.label), length(dropped)))

# ---- trait vector ---------------------------------------------------------
hybrid_counts <- setNames(sp$Descendants, sp$SpeciesName)
common <- intersect(tree_species$tip.label, names(hybrid_counts))
tree_phylo <- drop.tip(tree_species, setdiff(tree_species$tip.label, common))
hybrid_counts <- hybrid_counts[common]

# ---- branch lengths -------------------------------------------------------
tree_phylo$edge.length <- rep(1, nrow(tree_phylo$edge))
tree_grafen <- compute.brlen(tree_phylo, method = "Grafen")

# ---- Table 1, species-level row -------------------------------------------
K_result      <- phylosig(tree_grafen, hybrid_counts, method = "K",      test = TRUE)
lambda_result <- phylosig(tree_grafen, hybrid_counts, method = "lambda", test = TRUE)

trait_df <- data.frame(hybrids = hybrid_counts[tree_grafen$tip.label])
rownames(trait_df) <- tree_grafen$tip.label
tree4d <- phylo4d(tree_grafen, trait_df)
ps <- phyloSignal(tree4d, methods = c("K", "I", "Cmean"), reps = 999)

cat("\n--- Table 1, species-level row ---\n")
cat(sprintf("Blomberg's K       %.3f  (p = %.3f)\n", K_result$K, K_result$P))
cat(sprintf("Pagel's lambda     %.3f  (p = %.3g)\n",
            lambda_result$lambda, lambda_result$P))
cat(sprintf("Moran's I          %.3f  (p = %.3f)\n", ps$stat$I,     ps$pvalue$I))
cat(sprintf("Abouheif's Cmean   %.3f  (p = %.3f)\n", ps$stat$Cmean, ps$pvalue$Cmean))

species_row <- data.frame(
  tree = "species level",
  statistic = c("Blomberg K", "Pagel lambda", "Moran I", "Abouheif Cmean"),
  value = c(K_result$K, lambda_result$lambda, ps$stat$I, ps$stat$Cmean),
  p     = c(K_result$P,  lambda_result$P,     ps$pvalue$I, ps$pvalue$Cmean))

# ---- Figure 3 -------------------------------------------------------------
tree_clean <- tree_phylo
tree_clean$edge.length <- rep(1, nrow(tree_clean$edge))
tree_clean <- compute.brlen(tree_clean, method = "Grafen")

abbrev <- setNames(
  gsub(" ", "_", sub("^([A-Za-z])[a-z]+\\s+", "\\1. ", sp$Binomial)),
  sp$SpeciesName)
tree_clean$tip.label <- abbrev[tree_clean$tip.label]
counts_clean <- setNames(hybrid_counts, abbrev[names(hybrid_counts)])
counts_clean <- counts_clean[tree_clean$tip.label]

stopifnot(all(tree_clean$tip.label == names(counts_clean)))
stopifnot(!any(is.na(counts_clean)))

# Tree drawn from ape's layout so it can share the house ggplot style: a
# rectangular phylogram on the left, descendant counts as log-scale bars on the
# right, each subgeneric rank shaded as a band across both panels.
x_node <- node.depth.edgelength(tree_clean)
y_node <- node.height(tree_clean)
n <- Ntip(tree_clean)
ed <- tree_clean$edge
horiz <- data.frame(x = x_node[ed[, 1]], xend = x_node[ed[, 2]],
                    y = y_node[ed[, 2]], yend = y_node[ed[, 2]])
vert <- do.call(rbind, lapply(unique(ed[, 1]), function(p) {
  ys <- y_node[ed[ed[, 1] == p, 2]]
  data.frame(x = x_node[p], xend = x_node[p], y = min(ys), yend = max(ys))
}))

rank_of <- setNames(sp$LastSubgenus, abbrev[sp$SpeciesName])
tips <- data.frame(label = tree_clean$tip.label, x = x_node[1:n], y = y_node[1:n],
                   count = as.numeric(counts_clean[tree_clean$tip.label]))
tips$rank <- rank_of[tips$label]
tips$text <- sub("_", " ", tips$label)

ranks <- do.call(rbind, lapply(split(tips, tips$rank), function(d) {
  data.frame(rank = d$rank[1], ymin = min(d$y) - 0.5, ymax = max(d$y) + 0.5,
             ly = mean(d$y))
}))
ranks <- ranks[order(ranks$ymin), ]
ranks$fill <- ifelse(ranks$rank == "Cattleya", seq_blue[["light"]],
                     ifelse(seq_len(nrow(ranks)) %% 2 == 0, ink[["band"]], NA))

y_lim <- c(0.5, n + 0.5)
depth <- max(x_node)
bands <- function() geom_rect(data = subset(ranks, !is.na(fill)),
                              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
                                  fill = I(fill)), inherit.aes = FALSE)

p_tree <- ggplot() +
  bands() +
  geom_segment(data = rbind(horiz, vert), aes(x, y, xend = xend, yend = yend),
               colour = ink[["secondary"]], linewidth = 0.35, lineend = "square") +
  geom_text(data = tips, aes(depth + depth * 0.02, y, label = text),
            hjust = 0, size = 2.5, fontface = "italic", colour = ink[["primary"]]) +
  geom_text(data = ranks, aes(depth * 1.5, ly, label = rank), hjust = 0, size = 2.7,
            colour = ink[["secondary"]], fontface = "bold") +
  scale_x_continuous(limits = c(0, depth * 1.9), expand = expansion(0)) +
  scale_y_continuous(limits = y_lim, expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_cattleya() +
  theme(axis.text = element_blank(), panel.grid.major = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = margin(0, 0, 0, 0))

top <- head(tips[order(-tips$count), ], 6)
p_bars <- ggplot(tips) +
  bands() +
  geom_segment(aes(x = 0, xend = log10(count), y = y, yend = y),
               colour = seq_blue[["dark"]], linewidth = 1.6, lineend = "butt") +
  geom_point(aes(log10(count), y), colour = seq_blue[["dark"]], size = 1.1) +
  geom_text(data = top, aes(log10(count), y, label = scales::comma(count)),
            hjust = -0.2, size = 2.5, colour = ink[["secondary"]]) +
  # log10 applied by hand so the shaded bands can run from -Inf
  scale_x_continuous(breaks = 0:4, labels = scales::comma(10^(0:4)),
                     expand = expansion(mult = c(0.02, 0.2))) +
  scale_y_continuous(limits = y_lim, expand = expansion(0)) +
  labs(x = "Descendant hybrids (log scale)", y = NULL) +
  theme_cattleya() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = ink[["grid"]], linewidth = 0.35),
        plot.margin = margin(0, 4, 0, 6))

fig3 <- (p_tree | p_bars) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    title = "Heavily used parents cluster in the labiata group",
    subtitle = "Cattleya species on the subgeneric tree; the labiata group (rank Cattleya) is shaded blue",
    theme = theme_cattleya())
save_fig(fig3, file.path(out_dir, "fig3_hybrid_counts_on_species_tree.png"),
         width = 8.5, height = 12)

cat(sprintf("\nwrote %s/fig3_hybrid_counts_on_species_tree.png\n", out_dir))

# ---- Table 1, subgeneric-rank row -----------------------------------------
# Hybrid counts summed per terminal subgeneric rank, on the rank tree itself.
# Stellata is excluded as in the submitted analysis; Crispae is not a tip of
# the tree and drops out at the intersect. Run last so the seeded permutation
# p-values of the species-level row above are unchanged.
rank_counts <- tapply(sp$Descendants, sp$LastSubgenus, sum, na.rm = TRUE)
rank_counts <- rank_counts[names(rank_counts) != "Stellata"]
common_ranks <- intersect(tree_lowestrank$tip.label, names(rank_counts))
tree_rank <- drop.tip(tree_lowestrank,
                      setdiff(tree_lowestrank$tip.label, common_ranks))
tree_rank <- compute.brlen(tree_rank)          # Grafen, the ape default
rank_counts <- rank_counts[tree_rank$tip.label]

K_rank      <- phylosig(tree_rank, rank_counts, method = "K",      test = TRUE)
lambda_rank <- phylosig(tree_rank, rank_counts, method = "lambda", test = TRUE)
rank_df <- data.frame(hybrids = as.numeric(rank_counts))
rownames(rank_df) <- tree_rank$tip.label
ps_rank <- phyloSignal(phylo4d(tree_rank, rank_df),
                       methods = c("I", "Cmean"), reps = 999)

cat(sprintf("\n--- Table 1, subgeneric-rank row (%d ranks) ---\n",
            length(tree_rank$tip.label)))
cat(sprintf("Blomberg's K       %.3f  (p = %.3f)\n", K_rank$K, K_rank$P))
cat(sprintf("Pagel's lambda     %.5f  (p = %.3g)\n",
            lambda_rank$lambda, lambda_rank$P))
cat(sprintf("Moran's I          %.3f  (p = %.3f)\n", ps_rank$stat$I, ps_rank$pvalue$I))
cat(sprintf("Abouheif's Cmean   %.3f  (p = %.3f)\n",
            ps_rank$stat$Cmean, ps_rank$pvalue$Cmean))

rank_row <- data.frame(
  tree = "subgeneric rank",
  statistic = c("Blomberg K", "Pagel lambda", "Moran I", "Abouheif Cmean"),
  value = c(K_rank$K, lambda_rank$lambda, ps_rank$stat$I, ps_rank$stat$Cmean),
  p     = c(K_rank$P, lambda_rank$P, ps_rank$pvalue$I, ps_rank$pvalue$Cmean))
write.csv(rbind(rank_row, species_row),
          file.path(out_dir, "table1_signal_tests.csv"), row.names = FALSE)
cat(sprintf("wrote %s/table1_signal_tests.csv\n", out_dir))
