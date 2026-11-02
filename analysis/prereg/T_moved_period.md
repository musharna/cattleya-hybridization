# Pre-registration T: why moved species have more primaries (frozen 2026-10-27, before the moved-vs-rest split by period is seen)

## Why
In test C, species moved into *Cattleya* from other genera have about 2.5 times the primary hybrids per year since description (molecular tree, rate ratio 2.5, 1.5–4.2). Two readings:
- **Naming:** before the moves, crossing one of these species with a *Cattleya* made an intergeneric hybrid, and breeders sought those out. Then the excess should shrink after the move.
- **Appeal:** breeders like the species themselves. Then the excess should stay.

Seen before freezing: pooled counts over all species (pre 1,823, post 631 species-primary counts) and the combination years. Not seen: any count split by moved vs rest and period.

## Data (`analysis/derive_period_counts.py` → `data/species_period_counts.csv`)
- Primaries counted as in C (T1). Grexes without a registration year are dropped (13).
- **Moved:** WCVP v15 Cattleya combination published 2006–2014 (52 species). The 5 species moved before 1890, before registration began, count with the rest.
- **Cutoff:** each moved species' combination year (mostly 2008); 2008 for the rest. Pre = registered in or before the cutoff year, post = after.
- **Exposure:** pre = cutoff − description year, post = 2025 − cutoff (the same 2025 end as C). A row with exposure under 1 year is dropped and the number reported.

## Model
brms Poisson, one row per species × period, on the combined molecular tree (91 species):
`n ~ moved * post + offset(log exposure) + (1 | gr(species, cov = A)) + (1 | species_iid) + (1 | obs)`.
4 chains × 4,000 iterations, seed 1, Rhat < 1.01 required.
**Readout:** ratio = exp(moved:post) = rate ratio after / rate ratio before, with its 95% CrI.

## Verdict rules
- **Naming:** ratio CrI entirely below 1, and the placebo ratio's CrI includes 1.
- **Appeal:** ratio CrI includes 1 and lies inside [0.5, 2].
- **Can't tell:** anything else.
- **Placebo:** the same model on data up to the cutoff only, split at 1985 instead.
- **Sensitivity (reported, not decisive):** cutoff 2014 for all; and a post-period covariate log(1 + pre count) for saturation (each species pair can be registered only once).

## Planted truth first (gate)
Same tree, exposures and groups. Counts are simulated as Poisson with log-rate = period intercept + log exposure + moved effect + species effect + observation noise:
- Period intercepts are set so the expected pooled totals match 1,823 and 631.
- The species effect is √0.3·BM + √0.7·iid, scaled to sd(log1p T1) as in C; the observation noise has sd 0.3.
- **Constant:** moved effect log 2.5 in both periods (correct verdict: appeal).
- **Naming:** log 2.5 before, 0 after (true ratio 0.4; correct verdict: naming).
- 10 reps each through the production model.

**The test is usable only if ≥ 8/10 of each scenario get the correct verdict and ≤ 2/10 of the constant scenario are called naming.** If not, stop: the result is "can't tell", and the real split is not opened.

Outcome goes into `docs/v2_deeper.md` only after the user reviews it.
