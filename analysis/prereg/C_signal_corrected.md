# Pre-registration C: effort-controlled phylogenetic signal (frozen 2026-10-25, before trait values are computed)

## Why
v1 (`analysis/species_tree_signal.R:67,78`) tested `Descendants`, which is every hybrid with the species anywhere in its ancestry. That totals 443,474 over 100 species against 18,739 hybrids, and includes intergeneric hybrids. Pagel's λ was fit to the raw counts under a Gaussian model. Cumulative counts like this cluster on a tree almost by construction.

## Tree
Identical to v1: subgeneric topology with species polytomies, unit edges, Grafen. The tip set must equal v1's (asserted).

## Traits (from the 2025 scrape via `analysis/register_qc.py` resolution; non-synonym records only; a "var./subsp./f." form maps to its binomial)
- **T1 (primary):** registered primary grexes (both parents species) in which species s is a parent. A grex counts once per distinct parent species; selfings and two-forms-of-one-species count once.
- **T2 (direct use):** grexes of any type with s as an immediate parent.
- **Exposure E:** 2025 − description year from the species list (original description year, verified for C. coccinea 1836). Floor 1.
- **Transferred:** WCVP v15 basionym genus is not Cattleya (ex Laelia, Sophronitis, …). A missing basionym counts as not transferred, and the number of such cases is reported.

## Positive control (pipeline)
Reproduce v1 on `Descendants`: λ within 0.05 of 0.337, p < 1e-4.

## Tests
1. λ and K on log1p(T1), and on log1p(T2) as a secondary.
2. **Placebo:** λ on log(E), the description era, which no hybridization process sets.
3. **Effort-controlled:** brms Poisson on T1, with offset log(E), fixed effect `transferred`, a phylogenetic random effect (vcv of the Grafen tree, correlation) and an observation-level random effect. H2 = σ²_phylo / (σ²_phylo + σ²_obs). Report the posterior median and 95% CrI. 4 chains × 2,000 iterations, seed 1. Rhat < 1.01 is required, otherwise the result is reported as unconverged.
4. **Planted truth:** simulate log-rate = log(E) + √λ·BM + √(1−λ)·iid on the same tree (variance matched to the observed log1p(T1)), with Poisson counts.
   - 200 reps × λ ∈ {0, 0.3, 0.6}, through test 1.
   - 10 reps × λ ∈ {0, 0.6}, through test 3.
   - Report the false-positive rate at λ=0 (p<0.05) and the power at 0.3/0.6.

## Verdict rules (frozen)
- **Test invalid on this tree** if the FPR at planted λ=0 exceeds 0.10 for test 1. Test 1 then cannot carry the claim, and test 3 decides.
- **v1 claim NOT supported** ("clustering of hybridization") if any of:
  - (a) test 1 on log1p(T1) gives p ≥ 0.05 (when the test is valid);
  - (b) the test 3 H2 95% CrI lower bound is < 0.05;
  - (c) placebo λ ≥ observed λ on log1p(T1).
- **Supported** if none of (a)–(c) holds.
- Either way, the `transferred` coefficient and its CrI are reported. If H2 drops by more than half when `transferred` is removed vs included, the signal is called partly nomenclatural.
- Outcome goes to `CORRECTIONS.md` only after the user reviews it (Q4).

## Amendment 1 (2026-10-25, after the planted pilots only; no real-data result seen)
- **Test 1 simulation (200 reps):** the false-positive rate at planted λ=0 was 1.00 (mean λ̂ 0.36), because the exposure term is itself structured on the tree. Under the rule above, test 1 is invalid and test 3 decides.
- **Iterations raised to 4,000 (2,000 warmup).** In a one-rep brms pilot the planted λ=0 fit reached Rhat 1.017, with low ESS and 1 divergence. That breaks the Rhat < 1.01 rule at 2,000 iterations.
- **Test 3 simulations run as 4 parallel shards** (λ ∈ {0, 0.6} × 2 shards × 5 reps = 10 per λ), with distinct seeds 1000·shard + 10λ. Same model, same totals.
- **Wording fix to the "partly nomenclatural" rule:** it fires if H2 with `transferred` included is less than half of H2 without it, i.e. including `transferred` absorbs more than half the phylogenetic signal. The rule as first written had the comparison backwards.

## Amendment 2 (2026-10-25): re-run on a molecular tree. Written after the real result on the v1 tree, before any molecular tree exists
**Why:** on the v1 tree (10 subgeneric nodes, species as polytomies, unit edges) test 3 gave H2 0.35 [0.005, 0.85], meaning it could not tell. Its power came from planted data on that same weak tree. A resolved tree may decide. The trait values T1, E and `transferred` are already known and do not change. Only the tree changes, and it is built from sequences alone, without reference to the traits.

**Data (GenBank, fetched once, accession list saved in `data/`):** ITS (nuclear) and matK and trnL-F (plastid) records under Cattleya and its former genera (Laelia, Sophronitis, Hadrolaelia, Hoffmannseggella and similar), mapped to WCVP v15 accepted names. Probe 2026-10-25: ITS covers 92 of the 100 species, a plastid marker covers about 45, and 55 have both.
- One accession per species per marker, the longest. Minimum length: ITS 400 bp, matK and trnL-F 500 bp. Records named as hybrids (×) are excluded.
- Outgroups: Laeliinae genera outside Cattleya with ITS (Guarianthe, Brassavola, Rhyncholaelia; ≥ 2 of them). They root the tree and are then pruned.

**Trees:**
- **Primary: combined.** ITS + matK + trnL-F concatenated, with plastid missing where absent. mafft `--auto` per marker. IQ-TREE partitioned by marker, ModelFinder, 1,000 ultrafast bootstraps, seed 1.
- **Sensitivity:** ITS-only and plastid-only trees, built the same way.
- **Ultrametric:** `ape::chronos` (λ = 1) on the rooted ML tree, scaled to height 1, then the vcv used as in test 3.
- **Species set:** the primary analysis uses the species in the combined tree (expected 92). Sensitivity: all 100, with the 8 unsequenced species grafted as polytomies at the crown of their lowest rank in the v1 classification, edge length 0.
- **Sanity check, not a gate:** recovery of Van den Berg (2014) subgenera as clades. This is not independent, because that classification was built from trees on these same markers.

**Tests and rules:** tests 1–4 and the verdict rules (a)–(c) as frozen above, with the Amendment 1 settings (4,000 iterations, Rhat < 1.01), run on the primary tree.
- The test 1 validity simulation (200 reps × λ ∈ {0, 0.3, 0.6}) and the test 3 planted truth (10 reps × λ ∈ {0, 0.6}) are re-run on the new tree.
- Test 3 is called **invalid** on the new tree if more than 2 of 10 planted λ=0 fits pass the bar. It is called **underpowered** if fewer than 8 of 10 planted λ=0.6 fits pass.
- **Tree uncertainty:** test 3 is also fit on 20 bootstrap trees (the first 20 in IQ-TREE's `.ufboot` file, each ultrametricized the same way), with posterior draws pooled. The pooled H2 is reported next to the ML-tree H2.
- **"Tree-dependent":** if the verdict on the primary ML tree differs from the pooled-bootstrap verdict, or from either single-marker tree, the result is labelled tree-dependent and not supported as a robust claim.
- The v1-tree result stays reported alongside. `CORRECTIONS.md` changes only after the user reviews it.

## Amendment 3 (2026-10-25): rooting. Trees built, no H2 or other trait result on them seen yet
Amendment 2 assumed the outgroups would form a clade. They don't.
- **Plastid tree:** the outgroups form a clade.
- **Combined tree:** one tip breaks it. *Brassavola acaulis* (ITS only, AY008590) sits inside a *Cattleya* clade at 78% UFBoot; without it, *Cattleya* is a clade.
- **ITS-only tree:** removing any one or two outgroup tips does not make *Cattleya* a clade.
- Two GenBank fixes were made before this: 6 rpl32–trnL(UAG) records had been taken as trnL-F, and records whose title does not name the marker are now excluded.

The rule for every tree:
1. If the study species form a clade, root on all outgroups.
2. Otherwise, drop the smallest set of at most 2 outgroup tips that makes them a clade (combined: *B. acaulis*), list the dropped tips, and root on the remaining outgroups.
3. If no such set exists (ITS-only), root on *Rhyncholaelia*, the outgroup genus recovered as a clade in that tree, and record that ITS alone does not support *Cattleya* being a clade here.

The same rule applies to each bootstrap tree. A bootstrap tree that needs step 3 is rooted on *Rhyncholaelia* too, and the number of such trees is reported.

## Amendment 4 (2026-10-26): power at weak signal. Written after the real molecular result (H2 0.19, 0.004–0.49), before any λ 0.3 test 3 fit exists
**Why:** test 3 was planted only at λ 0 and 0.6. The real interval excludes strong signal but not weak signal, so whether weak signal is absent or merely undetectable depends on test 3's power at an intermediate value, which was never measured.

**Run:** combined molecular tree, `C_TREE=../data/molecular_combined.nwk C_TAG=mol Rscript signal_corrected.R sim3 5 0.3 SHARD` for shards 1 and 2 (10 reps, seeds 1000·shard + 3), model and settings unchanged.

**Readout:** a fit detects signal if its H2 95% lower bound is ≥ 0.05 (the bar of rule (b)). Fits with Rhat ≥ 1.01 are listed and left out of the denominator, which is reported.
- **≥ 8 of 10 detect:** weak signal at this level would usually have been seen, so the real result is read as little or no signal.
- **≤ 4 of 10 detect:** weak signal at this level is undetectable with these data; the write-up says it can be neither shown nor ruled out.
- **5–7:** reported as partial power, with the count.

With 10 reps the power estimate is coarse (8/10 has a 95% binomial interval of about 0.44–0.97); the count is reported, not a rate. The median recovered H2 at planted λ 0.3 is reported next to the real 0.19, descriptively only. `CORRECTIONS.md` and the README change only after the user reviews the result.

## Amendment 5 (2026-10-27): exposure bug. Found while building prereg T; written before any rerun
The species list's `Year` (used for E since the first freeze) is the year of the current *Cattleya* name, not the original description. For the 52 species moved in 2006–2014 it is 2006–2014, so their exposure was about 17 years instead of a century or more. (The one check in the first freeze, *C. coccinea* 1836, was described in *Cattleya* and so could not show this.) For 50 of the 52 the year differs from WCVP; for the other 48 species WCVP agrees exactly in 40, and the 8 differences are 1–36 years (described first as a variety or under an older name).
- **Fix:** E = 2025 − `orig_year`, the WCVP v15 publication year of the basionym (or of the name itself if it has none). *C. gloedeniana* has no WCVP year and keeps 2008.
- **Rerun, tag `mol_e5`:** test 3 on the combined molecular tree (real), and the test 3 planted truth at λ 0, 0.6 and 0.3 (exposure enters the simulation). Same model, rules and bars as before. The v1 tree is rerun too (tag `e5`).
- **Not rerun now:** bootstrap, ITS-only, plastid-only and 100-species fits, and test 1 (already invalid). They are rerun only if the combined-tree verdict changes.
- Every earlier number that uses E (all of test 3, the `transferred` rate ratio, the λ 0.3 power) is superseded by the rerun. That includes the public "2.5 times" claim.
- **Addendum (after the combined-tree rerun; verdict unchanged):** the bootstrap, ITS-only, plastid-only and 100-species fits are rerun anyway (tags `mol_e5`, `mol_its_e5`, `mol_plastid_e5`, `mol_all100_e5`). The rule above did not require it, but the published figure and text cite those numbers, and they used the faulty exposure. They are reported; the verdict rules are unchanged.

## Amendment 6 (2026-10-28): refit the unconverged sensitivity fits. Written after the e5 results, before any refit
Seven published e5 fits have Rhat ≥ 1.01 (with `transferred`): ITS-only 1.020, 100-species 1.016, and bootstrap trees 1, 4, 7, 15 and 17 (1.010–1.027). They are reported as "did not fully converge".
- **Refit, tag `e6`:** the same fits, model, data and seed, at 8,000 iterations (4,000 warmup), thin 2 (so each fit keeps 8,000 draws, as before) and adapt_delta 0.99. ITS-only and 100-species are rerun in `real` mode; the five bootstrap trees are refit one at a time.
- **Readout, fixed now:** each refit replaces its e5 value in the figure and text, whatever it shows. The pooled bootstrap uses the e6 draws for those five trees and the e5 draws for the other 15. A refit still at Rhat ≥ 1.01 keeps its label; there is no second attempt.
- Verdict rules unchanged, including Amendment 2's "tree-dependent" rule. No other fit is refit.

## Amendment 7 (2026-10-28): refit the unconverged planted-truth fits. Written after the e5 sims, before any refit
Six e5 planted-truth fits on the combined tree have Rhat ≥ 1.01: at λ 0, shard 1 reps 1 and 4 and shard 2 reps 1, 4 and 5 (1.012–1.018); at λ 0.3, shard 2 rep 3 (1.028). The λ 0.6 fits all converged.
- **Refit, tag `mol_e7`:** reps within a shard share one random-number stream, so one rep cannot be redrawn alone. The three shards that hold them (λ 0 shards 1 and 2, λ 0.3 shard 2) are rerun whole, with the same seeds and model, at the Amendment 6 sampler settings (8,000 iterations, 4,000 warmup, thin 2, adapt_delta 0.99).
- **Readout, fixed now:** all five reps of each rerun shard replace their e5 values, whatever they show. λ 0.3 shard 1 and both λ 0.6 shards keep e5. The rules of Amendments 2 and 4 are unchanged: a fit detects signal if its H2 95% lower bound is ≥ 0.05; a fit still at Rhat ≥ 1.01 is listed and left out of the denominator; there is no second attempt.
- **Same-data check:** rep 1 of each shard is drawn before any fit, so its data are the same by construction. If brms does not consume R's random stream, reps 2–5 are too. For the nine e5 reps in these shards that converged, the rerun H2 median should match within 0.05. If any differs by more, the rerun is reported as new draws, not a refit; the readout above still applies.
- `CORRECTIONS.md`, the README and `docs/v2_deeper.md` change only after the user reviews the result.

## Amendment 8 (2026-10-28): power curve. Written after the e7 planted-truth counts (λ 0: 0/10, 0.3: 5/10, 0.6: 8/10), before any new fit and before any curve is fit
**Why:** the public text says weak clustering "can neither be shown nor ruled out". A power curve turns that into a bound: the weakest planted λ the test detects 80% of the time (the bar of Amendment 2's "underpowered" rule).
- **New fits, tag `mol_e8`:** combined molecular tree, planted λ 0.45 and 0.8, 2 shards × 5 reps each (seeds 1000·shard + round(10λ), R rounding: 1004 and 1008, 2004 and 2008), at the Amendment 6 sampler settings. Model unchanged.
- **Data for the curve:** these 20 fits plus the 30 e5/e7 fits at λ 0, 0.3 and 0.6 as published after Amendment 7. A fit detects signal if its H2 95% lower bound is ≥ 0.05. Fits with Rhat ≥ 1.01 are listed and left out; there is no refit.
- **Curve:** logistic regression of detect on λ (glm, binomial, logit). λ80 = the λ where the fitted probability is 0.8. 95% interval: 2,000 bootstrap resamples of fits within each λ, percentile interval, seed 1; resamples where the slope is not positive or λ80 falls outside [0, 1] are counted and reported, and set to 1 (no 80% power within range).
- **Readout, fixed now:** report λ80 and its interval, the count per λ, and the median recovered H2 per λ (descriptive). If λ80 ≤ 0.8, the public sentence becomes: clustering at λ ≥ λ80 would have been detected about 4 times in 5; weaker clustering can be neither shown nor ruled out. If λ80 > 0.8, it says the test reaches 80% power only above λ 0.8. Public text changes only after the user reviews it.
- Verdict rules for the real data are unchanged.
