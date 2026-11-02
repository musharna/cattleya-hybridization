# Pre-registration W: do breeders cross close or distant species pairs? (2026-10-29)

Written before any pair-level count has been computed. Nothing below has been looked at except the per-species totals already published in C.

**Question.** Tests C, U and V each asked about single species, and each hit the same limit of 91 species. Here the unit is the species pair. Among registered primary hybrids, are closely related pairs crossed more or less often than distant ones, once each species' overall popularity and how long both species have been known are allowed for? Only successful crosses get registered, so the answer is about registered crossing: breeders' choice and cross success together. It cannot separate the two.

**Prior art** (map in project memory, `cattleya_expansion_prior_art_2026-10-29`). Experimental orchid crosses (Scopece 2007, 2008), natural-hybrid records (Brown 2023, Saward 2025) and Masters 2025 treat a cross as a readout of compatibility. None models breeders' choice of pairs over a register with exposure. The method is borrowed from trade economics: a Poisson pseudo-likelihood gravity model (Silva & Tenreyro 2006). Inference uses node-label permutation (QAP) from network sociology.

**Data.**
- **Species and pairs:** the 91 species on `data/molecular_combined.nwk`, giving 4,095 unordered pairs.
- **Outcome:** `n_ij` is the number of registered primary grexes whose two parents are species i and j, under T1 rules as in C (`register_qc` resolution; synonyms dropped; intraspecific crosses excluded). Derived by `analysis/derive_pair_counts.py` into `data/species_pair_counts.csv`, which holds species names only.
- **Exposure:** `E_ij` = 2025 − max(orig_year_i, orig_year_j), with a minimum of 1. orig_year is the WCVP basionym year, as in C after Amendment 5.
- **Distance:** `d_ij` is the cophenetic distance on the tree, standardized to mean 0 and SD 1 over the pairs.

**Model.** Poisson GLM: `n_ij ~ β·d_ij + a_i + a_j + offset(log E_ij)`. Each pair gets the sum of its two species' effects, entered as a 0/1/1 indicator design with one species as reference. The species effects absorb each species' popularity. With Poisson MLE, fitted totals per species equal observed totals; the script asserts this as the mechanism check.
- **Estimate:** β, reported as RR = exp(β) per SD of distance.
- **95% CI:** Wald, using a dyadic-robust variance (pairs sharing a species are correlated).
- **p-value:** permute the species labels on the distance matrix (999 permutations), refit, and take the two-sided p as the share of |β_perm| ≥ |β|.

**Planted truth, run first, through the same script.**
- **Simulated counts:** `n_ij ~ Poisson(E_ij · exp(c + a_i + a_j + β·d_ij + ε_ij))`, with ε ~ N(0, 1) pair noise.
- **Species effects:** a_i = log((T1_i + 0.5)/E_i), centered. This uses only the published per-species totals.
- **Scaling:** c is set so that the expected total equals the number of species-species primaries among the 91 species. That number is printed by the derivation script before any pair is modelled.
- **Arms:** 100 simulations each at RR 0.67, at RR 1.5, and at RR 1 (null), with 199 permutations per simulation. Seed 3000 + replicate (null), 4000 + (RR 0.67), 5000 + (RR 1.5).
- **Overdispersion sensitivity:** the null and RR 0.67 arms repeated with ε SD 0.5 and SD 1.5 (50 simulations each). Reported, not gating.

**Gates (frozen):**
- **Power:** at least 80 of 100 runs detect the planted effect (p < 0.05, in the planted direction) in each planted arm. Otherwise the test is **underpowered** for that direction.
- **False positives:** at most 8 of 100 null runs give p < 0.05. Otherwise the test is **invalid**.
- **Coverage:** at least 85 of 100 planted 95% CIs contain the planted β. Otherwise the CI is reported as **uncalibrated**, and only the permutation p is read.

**Readout on the real data.** RR per SD of distance, with its 95% CI and permutation p.
- **closer pairs crossed more:** p < 0.05 and RR < 1.
- **more distant pairs crossed more:** p < 0.05 and RR > 1.
- **no strong preference:** p ≥ 0.05, CI within [0.67, 1.5], and neither direction underpowered.
- **unresolved:** anything else.

**Sensitivity (reported, not gating):**
- (a) pairs that have a WCVP nothospecies dropped, since some registered primaries remake natural hybrids (E gate);
- (b) the 20 bootstrap trees, giving the spread of RR;
- (c) the ITS-only tree.

**Secondary (reported, not gating).** The same model with counts split at 1960 (eras as in prereg V), adding a `d × late` term and per-era species effects. It asks whether the distance preference changed when the mini-catt breeding began.

**Second route.** The real β point estimate is refit in Python (statsmodels GLM, same design). It must agree with R to within 0.01.

**Figure.** Observed vs fitted crossing rate by distance decile, with the permutation null band.

`docs/v2_deeper.md` changes only after the user reviews the result.

## Amendment 1 (2026-10-29, before any simulation or real fit)
- **Species with zero registered primaries among the 4,095 pairs are dropped.** Their effect has no finite MLE (minus infinity), and their pairs carry no information about β once species effects are in. This is standard in Poisson gravity fitting. The real-data count of species kept is reported.
- **Faster fitting engine.** The model is fit by damped proportional fitting of the species effects plus Newton steps for β (`fit_fast`), because `glm.fit` takes 4 s per fit (100,000 fits ≈ 4.5 days). This is the same MLE: on three planted data sets it agrees with `glm.fit` to 1e-6 (`species_pairs.R selftest`).
- **Sensitivity-arm seeds.** Seed = base + replicate + 1000·round(10·SD), where base is 3000 (null) or 4000 (RR 0.67).

## Amendment 2 (2026-10-29)
- **Second route engine (decided before the real fit).** The second route uses a numpy IRLS instead of statsmodels, which fails to import against the installed scipy (`_lazywhere`). On a simulated design it agrees with R to 1e-4.
- **Grexes dated before a parent's description (decided after the main estimate was seen; affects only the secondary era split).** The first real run stopped at a check in the era step. One grex is dated 1902 with parents *C. liliputana* × *C. tenebrosa*, but *liliputana* was first described in 1973. The date is wrong, most likely a reused grex name carrying an older registration date. Rule, applied to every grex: a grex dated before either parent species' first description stays in `n`, but is left out of `n_early`/`n_late`. 13 grexes are affected. All 1,098 pair totals are unchanged, so the main estimate and sensitivities (a)–(c) are unchanged.
