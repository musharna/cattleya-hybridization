# Pre-registration X: did the mini-catt founder species spread between breeders along the breeder network? (2026-10-29)

Written before the breeder network has been built or joined to adoption. Only the outcome margins below have been looked at.

**Question.** After 1960 the mini-catt founder species (*C. coccinea*, *walkeriana*, *aclandiae*, *intermedia*) went from 1% to 16% of ancestry (analysis A). Breeders can be linked by who used whose grexes as parents. Did breeders start using these species directly sooner when breeders they were linked to had already done so?

**Prior art** (project memory, `cattleya_expansion_prior_art_2026-10-29`). Rose registrations have been counted by nursery (Moser & Rhode 2011). Dog-breed popularity behaves like random copying (Herzog 2004). Pedigree networks have been built with cultivars as the nodes (Pincot 2021). None follows a new parent species from breeder to breeder. The method is borrowed from network-based diffusion analysis (Franz & Nunn 2009), in a discrete-time event-history form.

**Data** (Originator field joined back to the March 2025 scrape by grex name). Breeder-level tables stay local in `~/orchid-data/cattleya_deeper/breeders/` because they name people; the repo gets aggregates only.
- **Breeder:** the Originator string, normalized to lowercase letters. Missing values (376 grexes) and "O/U" (567) are dropped. No fuzzy merging.
- **Adoption year** of breeder b: the first registration year of a grex originated by b that has a mini-catt founder species as a *direct* parent (T1 resolution rules). Adoption therefore cannot happen by inheritance through a parent grex.
- **Contacts** of b before year t: other breeders whose grexes b used as a parent, in grexes b registered before t. **Adopted contacts:** those who had adopted before t. Breeders who adopted before 1950 (66) count as adopted contacts.
- **Risk set:** breeder-years 1950–2024 in which b registered at least one grex and had not yet adopted. Breeders who adopted before 1950 are not at risk. Looked at so far: 379 breeders first used a founder species in 1950 or later, out of 1,958 breeders active from 1950.

**Model.** Logistic regression on breeder-years: `adopt ~ ns(year, 4) + log1p(grexes that year) + log1p(grexes before) + log1p(contacts) + β·log1p(adopted contacts)`. β is reported as an odds ratio per unit of log1p(adopted contacts), with a Wald 95% CI and p.

**Planted truth, run first, on the real network and risk set.**
- **Simulation:** adoption is simulated year by year. Exposure is recomputed from the simulated adoptions; pre-1950 adopters are kept as they are.
- **Baseline hazard:** taken from the asocial model, which is the model above without the β term, fit to the real data. That fit touches the outcome, but β is never fit before the gates are read.
- **Arms:** 100 simulations each at planted OR 1 (null), 1.5 and 2. Seeds 6000 + replicate, 7000 + replicate and 8000 + replicate.
- **Confound arms** (no social effect; adoption driven by a breeder trait shared along the network): trait z_b = standardized(own N(0,1) + mean of the breeder's contacts' noise). Adoption log-odds get + γ·z_b, with γ = 0.5 and γ = 1.0. 100 simulations each, seeds 9000 + replicate and 9500 + replicate.

**Gates (frozen):**
- **Power:** at least 80 of 100 runs at OR 1.5 give p < 0.05 with β > 0. Otherwise the test is **underpowered** at 1.5, and power at OR 2 is reported.
- **False positives:** at most 8 of 100 null runs give p < 0.05. Otherwise the test is **invalid**.
- **Confound (decides the wording, not validity):** if more than 8 of 100 runs flag at γ = 0.5, then a positive result cannot be called copying. It is worded "adoption clustered on the breeder network; copying and shared background cannot be told apart".

**Readout on the real data.** OR, 95% CI and p.
- **followed the network:** p < 0.05 and OR > 1. Called **consistent with copying** only if the γ = 0.5 confound arm flags at most 8/100 times.
- **no network effect detected:** p ≥ 0.05, the CI upper bound is below 1.5, and the test is not underpowered.
- **unresolved:** anything else.

**Description (not tested).**
- Yearly adopters with a Bass fit (innovation p, imitation q).
- The share of all post-1950 mini-catt grexes made by the first 10 post-1950 adopters.
- Per-species first-use curves.

**Sensitivity (reported, not gating):** (a) Registrant instead of Originator; (b) adoption of *C. walkeriana* alone; (c) exposure as any adopted contact (0/1).

`docs/v2_deeper.md` changes only after the user reviews the result.

## Amendment 1 (2026-10-29, after the gates; no real estimate had been printed)
- **Contact-only breeders.** The Registrant sensitivity stopped at a check: one contact appears only on undated or synonym grexes, so it was missing from the breeder list. Breeders now include every contact; such breeders are never at risk and never adopt. The Originator tables, which the main analysis and all simulations use, are byte-identical before and after.
- **Reused grex names.** 217 grex names (223 for Registrant) carry two different breeders, so those grexes get no breeder and are dropped, about 1% of the register.
