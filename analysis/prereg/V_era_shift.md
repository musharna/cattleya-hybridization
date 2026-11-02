# Pre-registration V: did the change in breeders' species use follow the phylogeny? (2026-10-29)

Written before any per-species count has been split by era. Analysis A already showed the shift in founder ancestry (the mini-catt founders rose from 1% to 16% of ancestry between the 1950s and 1980s, the labiata group fell from 67% to 37%). This asks a different question: is the change in each species' rate of primary hybrids, from before 1960 to after, clustered on the phylogeny?

**Data** (`analysis/derive_era_counts.py` → `data/species_era_counts.csv`, species-level only). The 91 species on the combined molecular tree. Primary grexes counted by the T1 rules (`derive_species_traits.py`) and split by registration year: **early** ≤ 1959, **late** ≥ 1960 (Morel's mericloning, and where A's founder mix widens). Exposure per era = years the species had been described within it: early E = 1960 − max(orig_year, 1850), late E = 2025 − max(orig_year, 1960), each at least 1. A species described in 1960 or later has a late row only. No other era split is run.

**Model.** brms Poisson, one row per species × era:
`T ~ late + transferred + offset(logE) + (1 + late || gr(species, cov = A)) + (1 + late || sp) + (1 | obs)`, where `sp` is the species as an independent factor. H2_shift = σ²(phylogenetic `late`) / (σ²(phylogenetic `late`) + σ²(independent `late`)). Same tree and seed as prereg C; sampler as C Amendment 6. Rhat < 1.01 required, otherwise the fit is reported as unconverged; the model is not respecified after the data are seen.

**Planted truth, run first, same script.** Intercepts simulated as in C (`sim_counts`, planted λ 0.3). The species' late shift = 1.0 × (√λs · BM + √(1−λs) · noise), standardized as in C; the era effect = log(total late rate / total early rate) from the real totals only. λs = 0.6 (2 shards × 5) and λs = 0 (2 shards × 5); seeds 1000·shard + 70 and + 80. A fit passes if the H2_shift 95% lower bound is ≥ 0.05.
- **Invalid** if more than 2 of 10 converged λs = 0 fits pass. **Underpowered** if fewer than 8 of 10 converged λs = 0.6 fits pass.

**Mechanism check on the real data:** at least 3 of the 4 mini-catt founders of A (*C. coccinea*, *walkeriana*, *aclandiae*, *intermedia*) must have a positive posterior-median shift. If not, the era split is not measuring the shift A found, and the H2_shift result is reported as not interpretable.

**Readout.** H2_shift median and 95% CrI. **Supported** if the lower bound is ≥ 0.05 and the test is neither invalid nor underpowered; **not supported** otherwise, with the reason. Per-species shifts (posterior medians) are plotted in tree order, descriptively.

`docs/v2_deeper.md` changes only after the user reviews the result.

## Amendment 1 (2026-10-29): two clarifications. Written after the era table was built (totals only printed), before any fit
- **Input:** the 8 March 2025 scrape (`Desktop/Code/bioinf/orchidroots_cattleya_hybrids.csv`), the file that reproduces `species_traits.csv`: the era counts sum to 2,453 against T1's 2,478, no species exceeds its T1, and the 25 missing are accounted for by the 13 undated grexes. A 10 April 2025 copy in Downloads does not reproduce T1 (4 species exceed it) and is not used.
- **"Shift" in the mechanism check** means the species' deviation from the average change (its phylogenetic plus independent `late` random effects), not the raw change including the era fixed effect. A's result is about shares of ancestry, so the relative change is the matching quantity.
