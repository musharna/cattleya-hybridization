# Data

This directory holds derived tables only; the raw records are not distributed here.

## Files

| File                                | Rows   | What it is                                                                                                                |
| ----------------------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------- |
| `hybrids_by_year_type.csv`          | 474    | Registered hybrids counted by year and hybrid type. Reproduces Figure 1.                                                   |
| `ancestors_by_year.csv`             | 18,718 | `(Year, Ancestors)` pairs, one per registered hybrid. No names, parentage, registrants or links. Reproduces Figure 2.      |
| `hybrids_per_year.csv`              | 158    | Hybrids registered per year, plus the running cumulative total. Reproduces Figures S1 and S2.                              |
| `species_hybrid_counts.csv`         | 100    | Per-species hybrid descendant counts and terminal subgeneric rank. This is the trait the phylogenetic signal tests run on. |
| `cattleya_subgeneric.nwk`           |        | Subgeneric-rank topology, hand-built from Van den Berg (2014).                                                             |

All four derived tables are produced by [`../analysis/derive_tables.R`](../analysis/derive_tables.R),
which is the only script that reads the raw records.

## Why derived only

The underlying records come from [orchidroots.com](https://www.orchidroots.com/), which republishes
the Royal Horticultural Society's International Orchid Register. Republishing a large part of that
register is not needed to reproduce this analysis, so these tables carry no grex names, parentage,
registrants, originators or links. For the records themselves, use orchidroots.com or the RHS
register directly; both are updated continuously.

## Cleaning applied before derivation

Reproduced from the submitted analysis, in `derive_tables.R`:

- hybrids flagged "Synonym of" in `Parentage` are dropped;
- species whose `Binomial` contains `subsp.` or ` var.` are dropped;
- the year 2025 is excluded, being a partial year at scrape time;
- species with no subgeneric rank recorded are dropped from the trait table.

Hybrids are classified by the capitalisation of each parent's second word. Both parents species
gives Primary, both parents hybrids gives Highly Complex, and anything else gives Complex. After
cleaning this leaves 18,739 hybrids, of which 1,468 are primary, 5,005 complex and 12,266 highly
complex.
