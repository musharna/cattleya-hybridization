# Pre-registration U: what predicts breeders' use of a species (2026-10-29)

Written before any of the predictors below has been joined to the hybrid counts. The covariate table was built from WCVP alone and only its marginal counts were looked at (listed here).

**Question.** Test 3 of prereg C found that the number of primary hybrids a species has is not clustered on the phylogeny once time since description is allowed for. What else predicts it?

**Data.** The 91 species on the combined molecular tree (`data/molecular_combined.nwk`). Outcome T1, exposure E and `transferred` exactly as in prereg C after Amendment 5. Predictors from WCVP v15 (`analysis/derive_species_covariates.py` → `data/species_covariates.csv`), native and non-extinct records only:
- `brazil`: native to at least one Brazilian level-3 area (78 of 91).
- `wet`: WCVP climate "wet tropical" (37) vs "seasonally dry tropical" (54).
- `litho`: WCVP life form mentions lithophyte (16).
- `z_areas`: log number of native level-3 areas, standardized (57 species have one).

Left out: photo counts and awards (they follow use), flower size (no data source checked here has it), CITES (every *Cattleya* is Appendix II).

**Model.** brms Poisson: `T1 ~ transferred + brazil + wet + litho + z_areas + offset(logE) + (1|gr(species, cov=A)) + (1|obs)`, same tree, priors and seed as prereg C, sampler as C Amendment 6 (8,000 iterations, 4,000 warmup, thin 2, adapt_delta 0.99). Rhat < 1.01 required, otherwise the fit is reported as unconverged.

**Planted truth, run first, same script and model.** Counts simulated as in C (`sim_counts`, planted λ 0.3) plus log(2)·(x − mean x) for each predictor, i.e. rate ratio 2 per unit (per SD for `z_areas`): 2 shards × 5 reps. A null set with no predictor effect: 2 shards × 5 reps. Seeds 1000·shard + 50 (planted) and + 60 (null).
- A predictor is **recovered** in a planted fit if its 95% CrI lies above 1. If fewer than 8 of 10 converged planted fits recover it, it is labelled **underpowered** and a null result for it is not read as "no effect".
- A predictor is **flagged** in a null fit if its 95% CrI excludes 1. If more than 2 of 10 converged null fits flag it, the test is **invalid** for that predictor.

**Readout on the real data, per predictor:** rate ratio (posterior median) and 95% CrI.
- **matters:** CrI excludes 1 (and the predictor is valid).
- **no strong effect:** CrI within [0.5, 2] and the predictor is not underpowered.
- **unresolved:** anything else.

H2 with the predictors is reported next to C's 0.12, descriptively. The predictors are geographically and phylogenetically clustered (most non-Brazilian species are in the labiata group), so fixed effects and the phylogenetic term compete for the same variation; that is reported, not resolved.

`docs/v2_deeper.md` changes only after the user reviews the result.
