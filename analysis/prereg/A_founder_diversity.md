# Pre-registration A: founder-species composition and diversity of the register (frozen 2026-10-25, before any composition is computed)

## Pedigree
- Non-synonym records from the 2025 scrape. Parents are resolved by `register_qc.resolve`.
- A species parent (listed or unlisted; var./subsp. forms map to the binomial) is a founder.
- A grex parent passes on its own composition. A child is ½ seed + ½ pollen.
- An ambiguous, unresolved or missing parent contributes "unknown" mass, which is reported per cohort and not renormalised away.
- Records are processed in topological order; cycles = 0 per QC.

## Metrics, per registration-decade cohort (all grexes registered in that decade)
- Founder-species share p_s = mean expected contribution of species s over the cohort.
- f_e = 1 / Σ p_s² (effective number of founder species).
- f_a = effective number of ancestors (Boichard 1997, 10.1186/1297-9686-29-1-5), by iterative selection of the ancestor with the largest marginal contribution. Marginal = cohort weight reaching the ancestor without passing through an already-selected ancestor × (1 − the ancestor's genome fraction from already-selected ancestors). Stop when 99% of the known mass is explained.
- Ratio f_e / f_a (bottleneck signal).

## Gates (in order; a failure stops A)
1. **Second route:** OrchidRoots' own `ancestry_json` (274 Cattleya grexes in `~/orchid-data/orchidroots_grex_crawl.csv`, computed by OrchidRoots, not by us). For matched grexes, the L1 distance between our composition and theirs, on the union of species with our unknown mass excluded, must have median < 0.10. Report the share with L1 < 0.05.
2. **Planted truth through the same code:**
   - (a) toy pedigrees with exact answers: f_e = 2 for an all-X×Y register; f_a = 1 when the whole cohort descends from one grex;
   - (b) a simulated register (2,000 grexes over 100 years, 40 founder species, preferential reuse of parents) with a planted popular parent, which must rank #1 in marginal contribution;
   - (c) a planted shift: founder group G's share is raised from 5% to 30% after year 60, and the share trajectory must show it.
3. **Positive control, known history:** the cohort share of mini-catt founders {coccinea, walkeriana, aclandiae, intermedia} must be higher in the 1980s than in the 1950s, and the labiata-group share {labiata, mossiae, trianae, warscewiczii, dowiana} lower. If not, the pedigree walk is called broken.

## Tests (reported after the gates pass)
- **Kill A ("only restates species introduction"):** Spearman ρ between decade f_e and the cumulative number of founder species used by that decade is > 0.9.
- **Null for the bottleneck:** degree-preserving edge swaps within each decade (parent→child edges swapped only when both parents are strictly older than both children; 10 × edges swaps; 200 draws). The observed f_e/f_a is called a bottleneck beyond parent-use frequency only if it lies outside the null 95% band.
- **Shrug test:** "diversity declined" alone is not a finding. The reported finding is the per-decade founder composition shift with its shares, set against Lai 2024 (79 cultivars) and the rose shift in Leroy 2025.

## Interpretation limit
A parent grex is in practice one clone, so contributions are expectations, not transmitted genomes. The scrape carries no clone names, so clone-level reuse is invisible.

## Results on the 2025 scrape (2026-10-25)
- Gate 1 PASS: 238 grexes matched OrchidRoots' `ancestry_json`, median L1 0.001, 85% below 0.05. Shared layer: the register's parentage text, so this checks our parsing and walk, not the register.
- Gate 2 PASS: toy cases exact; planted popular parent ranked #1; planted shift recovered against a flat twin (`tests/test_founders.py`, 8 mutants killed).
- Gate 3 PASS: mini-catt share 0.012 (1950s) → 0.163 (1980s); labiata group 0.669 → 0.374.
- Kill not met: ρ(f_e, cumulative species) = 0.63 (0.37 for decades with ≥ 20 grexes).
- Bottleneck test: NOT supported. f_e/f_a lies in 1.00–1.02 in every decade and outside the null band in 1 of 16 decades (2000, at the edge). f_a ≈ f_e is structural here: with ≤ 109 founder species, species out-contribute any single grex.

## Amendment 1 (post-hoc, written after seeing the null; exploratory on this scrape)
Observed f_e lies below the edge-swap null band in every decade 1890–2020 (about half the null), and the largest single-grex contribution lies above it 1910–1970. Concentration came from which lineages were crossed with which, not from how often parents were reused. Because this was seen after the fact, it is **not** a result on the 2025 scrape. It is frozen here as the confirmatory test for the rerun on the OrchidRoots Research Data List: observed f_e below the 2.5% null quantile in ≥ 10 of the decades 1910–2010, with the same code, the same null (200 draws, seeds 1000–1199) and the same decade cohorts.
