# Pre-registration E gate: natural vs registered hybrids, go/no-go (frozen 2026-10-25, before any pair table is computed)

The gate decides whether the full pair-level model of E (P(natural hybrid) from crossability, range overlap, collection effort and rank, with crossed species effects) is worth building. It does not test that model.

## Data
- Natural hybrids: WCVP v15 `wcvp_names.csv`, genus Cattleya, `species_hybrid` = ×, rank Species, status Accepted (104 names), parents from `hybrid_formula`. Formula epithets map to accepted names through WCVP (synonym → `accepted_plant_name_id`). A formula that does not give exactly two parents that map to accepted Cattleya species is dropped and counted.
- Ranges: WCVP v15 `wcvp_distribution.csv`, native L3 areas only (introduced = 0, extinct = 0, location_doubtful = 0).
- Registered primaries: register records (2025 scrape) whose two parents both resolve to species via `register_qc.resolve`, mapped to WCVP accepted names. A primary counts as **independent** evidence of crossability only if its epithet is grex-style (capitalised) and differs, case-insensitively, from the epithet of any nothospecies with the same parent pair (RHS files remakes of natural hybrids under the nothospecies name).
- Pair universe: unordered pairs of register species that map to an accepted WCVP Cattleya species with at least one native L3 area.
- Rank strata: same lowest subgeneric rank (`LastSubgenus` in `data/species_hybrid_counts.csv`) vs different. Pairs with a species lacking a rank enter the pooled table only.

## Gate readouts (all three must pass for GO)
1. **Mechanism readout (cultivation removes premating barriers):** of nothospecies whose parents are both in the pair universe, the share with an independent registered primary is ≥ 75%.
2. **Trade origin:** of those nothospecies, the share whose parents share no native L3 area is ≤ 20%. Above 20%: NO-GO (the list is partly a trade artefact). 15–20%: GO with a flag.
3. **Power:** of crossable pairs (pairs with an independent registered primary), the share that share a native L3 area is ≤ 80%. Above 80%: NO-GO (no geographic contrast at L3).

## Reported with the gate (descriptive, no go/no-go weight)
- The 2×2 table nothospecies × independent registered primary, pooled and per rank stratum, with Fisher exact odds ratios and p-values.
- The allopatric share among WCVP "Artificial Hybrid" names, set against readout 2.
- Counts dropped at every mapping step.

## Planted truth (through the same functions, before the real run)
- Toy WCVP + toy register: a synonym in a formula maps to the accepted parent; reversed formula order gives the same pair; a remake registered under the nothospecies epithet is not independent, a differently named grex of the same pair is; var./subsp. parents map to the binomial.
- Toy ranges: a planted allopatric nothospecies is counted allopatric; introduced and doubtful records are ignored.
- Each readout is watched to flip on a fixture built to fail it.

## Interpretation limit
L3 areas are coarse (country or state level), so "sympatric" at L3 is not co-occurrence. The register says a cross was made and flowered, not that it was easy.
