#!/usr/bin/env Rscript
# Build the derived tables that the public analysis runs on.
#
# This script is the ONLY thing that touches the raw orchidroots scrape. It
# reproduces the cleaning and hybrid-classification logic from
# the submitted analysis, then writes aggregate tables that carry no grex
# names, parentage strings, registrants, or source links.
#
# Usage (requires the raw CSVs, which are not distributed with this repo):
#   Rscript derive_tables.R <dir-with-raw-csvs> <output-dir>

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "../data-raw"
out_dir <- if (length(args) >= 2) args[2] else "../data"

hybrids <- read.csv(file.path(raw_dir, "orchidroots_cattleya_hybrids.csv"))
species <- read.csv(file.path(raw_dir, "orchidroots_cattleya_species.csv"))

# --- Cleaning, verbatim from the submitted analysis -------------------------
species_filtered <- species[grepl("subsp\\. |\\svar\\.", species$Binomial) == FALSE, ]
hybrids_filtered <- hybrids[!grepl("Synonym of", hybrids$Parentage), ]
hybrids_filtered <- hybrids_filtered[hybrids_filtered$Year != 2025, ]

# --- Hybrid classification, verbatim from the submitted analysis ------------
classify_hybrid <- function(parentage) {
  parents <- unlist(strsplit(as.character(parentage), " × "))
  is_hybrid <- function(name) {
    words <- unlist(strsplit(name, " "))
    if (length(words) < 2) return(FALSE)
    return(grepl("^[A-Z]", words[2]))
  }
  if (all(sapply(parents, function(x) !is_hybrid(x)))) {
    return("Primary Hybrid")
  } else if (all(sapply(parents, function(x) is_hybrid(x)))) {
    return("Highly Complex Hybrid")
  } else {
    return("Complex Hybrid")
  }
}
hybrids_filtered$HybridType <- sapply(hybrids_filtered$Parentage, classify_hybrid)

# --- Derived table 1: hybrid counts by year and type (Figure 1) -------------
t1 <- as.data.frame(table(hybrids_filtered$Year, hybrids_filtered$HybridType))
colnames(t1) <- c("Year", "HybridType", "Count")
t1$Year <- as.numeric(as.character(t1$Year))
t1 <- t1[order(t1$Year, t1$HybridType), ]
write.csv(t1, file.path(out_dir, "hybrids_by_year_type.csv"), row.names = FALSE)

# --- Derived table 2: ancestor count per hybrid, year only (Figure 2) -------
# Two numeric columns. No names, no parentage, no registrant, no links.
t2 <- hybrids_filtered[!is.na(hybrids_filtered$Year) & !is.na(hybrids_filtered$Ancestors),
                       c("Year", "Ancestors")]
t2 <- t2[order(t2$Year, t2$Ancestors), ]
write.csv(t2, file.path(out_dir, "ancestors_by_year.csv"), row.names = FALSE)

# --- Derived table 3: hybrids per year, raw and cumulative (Figs S1, S2) ----
per_year <- as.data.frame(table(hybrids_filtered$Year))
colnames(per_year) <- c("Year", "Hybrids")
per_year$Year <- as.numeric(as.character(per_year$Year))
per_year <- per_year[order(per_year$Year), ]
per_year$Cumulative <- cumsum(per_year$Hybrids)
write.csv(per_year, file.path(out_dir, "hybrids_per_year.csv"), row.names = FALSE)

# --- Derived table 4: per-species hybrid descendants (Figure 3, Table 1) ----
# `X..Descendants` is the "# Descendants" column: how many registered hybrids
# carry that species in their ancestry. This is the trait the phylogenetic
# signal tests are run on.
sp <- species_filtered[!is.na(species_filtered$Subgeneric.Ranks) &
                       !is.na(species_filtered$`X..Descendants`), ]
sp$LastSubgenus <- unlist(sapply(
  strsplit(sp$Subgeneric.Ranks, "\n"),
  function(x) {
    x <- x[x != ""]
    if (length(x) == 0) return(NA_character_)
    tail(x, 1)
  }
))
sp <- sp[!is.na(sp$LastSubgenus) & sp$LastSubgenus != "", ]
t4 <- data.frame(Binomial = sp$Binomial,
                 LastSubgenus = sp$LastSubgenus,
                 Descendants = sp$`X..Descendants`)
t4 <- t4[order(-t4$Descendants), ]
write.csv(t4, file.path(out_dir, "species_hybrid_counts.csv"), row.names = FALSE)

cat(sprintf("hybrids_by_year_type.csv   %d rows\n", nrow(t1)))
cat(sprintf("ancestors_by_year.csv      %d rows\n", nrow(t2)))
cat(sprintf("hybrids_per_year.csv       %d rows\n", nrow(per_year)))
cat(sprintf("species_hybrid_counts.csv  %d rows\n", nrow(t4)))
cat(sprintf("\nhybrids after cleaning:    %d (from %d raw)\n",
            nrow(hybrids_filtered), nrow(hybrids)))
cat(sprintf("species after cleaning:    %d (from %d raw)\n",
            nrow(sp), nrow(species)))
cat("\nHybrid type totals:\n")
print(table(hybrids_filtered$HybridType))
