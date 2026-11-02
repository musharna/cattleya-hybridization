#!/usr/bin/env Rscript
# Ultrametric study-species trees from the IQ-TREE output (prereg C, Amendment 2).
#
# What: roots each ML tree on the Laeliinae outgroups, prunes them and every tip that is not a study
#   species, makes it ultrametric with ape::chronos (lambda = 1) and scales it to height 1. Does the
#   same for the first 20 bootstrap trees of the combined run, and writes a 100-species version with
#   the unsequenced species grafted at the crown of their lowest rank. Reports subgenus recovery.
# Usage, from analysis/: Rscript molecular_chronos.R TREE_DIR
# Output: ../data/molecular_{combined,ITS,plastid,combined_all100}.nwk, ../data/molecular_boot20.nwk,
#   ../data/molecular_trees.json (tips, outgroup monophyly, subgenus recovery).
suppressPackageStartupMessages({ library(ape); library(jsonlite) })
tdir <- commandArgs(TRUE)[1]
sp <- read.csv("../data/species_traits.csv", stringsAsFactors = FALSE)
sp$tip <- gsub(" ", "_", sp$Binomial)
og_re <- "^(Guarianthe|Brassavola|Rhyncholaelia)_"
set.seed(1)

# rooting rule, prereg C Amendment 3: all outgroups if the study species form a clade; else drop the
# smallest set of <= 2 outgroup tips that makes them one; else root on Rhyncholaelia
ingroup_clade <- function(t, og) {
  t1 <- root(t, outgroup = og[1], resolve.root = TRUE)
  is.monophyletic(t1, setdiff(t1$tip.label, og))
}
ultra <- function(t) {
  og <- grep(og_re, t$tip.label, value = TRUE)
  stopifnot(length(og) >= 2)
  dropped <- character(0); step <- 1L
  if (!ingroup_clade(t, og)) {
    cand <- c(lapply(og, identity), combn(og, 2, simplify = FALSE))
    fix <- Filter(function(d) ingroup_clade(drop.tip(t, d), setdiff(og, d)), cand)
    if (length(fix)) {
      dropped <- fix[[1]]; step <- 2L  # singles come first, so this is a smallest set
      t <- drop.tip(t, dropped); og <- setdiff(og, dropped)
    } else {
      step <- 3L
      og <- grep("^Rhyncholaelia_", og, value = TRUE)
      stopifnot(length(og) >= 1)
      t1 <- root(t, og[1], resolve.root = TRUE)
      if (!is.monophyletic(t1, og)) og <- og[1]
    }
  }
  mono <- step == 1L
  t <- root(t, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
  t <- drop.tip(t, setdiff(t$tip.label, sp$tip))
  t$edge.length[t$edge.length <= 0] <- 1e-6
  ct <- suppressWarnings(chronos(t, lambda = 1, quiet = TRUE))
  class(ct) <- "phylo"
  ct$edge.length <- ct$edge.length / max(node.depth.edgelength(ct))
  stopifnot(is.ultrametric(ct, tol = 1e-6))
  attr(ct, "og_monophyletic") <- mono
  attr(ct, "root_step") <- step
  attr(ct, "dropped") <- dropped
  ct
}

rank_recovery <- function(t) {
  r <- split(sp$tip[sp$tip %in% t$tip.label], sp$LastSubgenus[sp$tip %in% t$tip.label])
  r <- r[lengths(r) >= 2]
  sapply(r, function(s) is.monophyletic(t, s))
}

out <- list()
trees <- list()
for (k in c("combined", "ITS", "plastid")) {
  t <- ultra(read.tree(file.path(tdir, paste0(k, ".treefile"))))
  trees[[k]] <- t
  write.tree(t, sprintf("../data/molecular_%s.nwk", k))
  rr <- rank_recovery(t)
  out[[k]] <- list(tips = Ntip(t), root_step = attr(t, "root_step"), dropped_outgroup_tips = attr(t, "dropped"),
                   ranks_tested = length(rr), ranks_monophyletic = sum(rr), ranks = as.list(rr))
  cat(sprintf("%s: %d study tips; rooting step %d (dropped: %s); ranks monophyletic %d/%d\n",
              k, Ntip(t), attr(t, "root_step"), paste(attr(t, "dropped"), collapse = ","), sum(rr), length(rr)))
}

boot <- read.tree(file.path(tdir, "combined.ufboot"))[1:20]
bt <- lapply(boot, ultra)
class(bt) <- "multiPhylo"
write.tree(bt, "../data/molecular_boot20.nwk")
out$boot20_tips <- unique(sapply(bt, Ntip))
out$boot20_root_steps <- as.list(table(sapply(bt, attr, "root_step")))
cat("bootstrap rooting steps:", paste(names(out$boot20_root_steps), unlist(out$boot20_root_steps), sep = ":", collapse = " "), "\n")

# 100-species sensitivity: graft each unsequenced species at the crown of its lowest rank
t <- trees$combined
miss <- sp[!sp$tip %in% t$tip.label, ]
for (i in seq_len(nrow(miss))) {
  mates <- intersect(sp$tip[sp$LastSubgenus == miss$LastSubgenus[i]], t$tip.label)
  node <- if (length(mates) >= 2) getMRCA(t, mates)
          else if (length(mates) == 1) t$edge[t$edge[, 2] == which(t$tip.label == mates), 1]
          else Ntip(t) + 1
  h <- max(node.depth.edgelength(t)) - node.depth.edgelength(t)[node]
  t <- phytools::bind.tip(t, miss$tip[i], edge.length = h, where = node)
}
stopifnot(Ntip(t) == nrow(sp), is.ultrametric(t, tol = 1e-6))
write.tree(t, "../data/molecular_combined_all100.nwk")
out$grafted <- miss$tip
write_json(out, "../data/molecular_trees.json", auto_unbox = TRUE, pretty = TRUE)
