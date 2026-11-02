# Pre-registration B and D: simulate first, open real numbers only if the simulation shows power (frozen 2026-10-25, before either outcome is computed on the register)

The simulations may use register **structure** that earlier steps already exposed: child registration years, the number of dated grex parents per child, the species pair of each primary, and per-species primary counts (C's `primary_T1`). They do not read the outcomes: parent ages (B) and later use as a parent (D).

## B: mericloning (Morel 1960) and the generation interval
- **Outcome:** for each grex with at least one grex parent of known year, the mean of (child year − parent year) over its grex parents.
- **Statistic:** child-level hinge regression y = β0 + β1·t + β2·(t − b)₊. F(b) = (SSE_linear − SSE_hinge) / (SSE_hinge / (n − 3)), with b on 1930–1990 in 1-year steps.
- **Pass (the break is real):** F(1960) > max F(b) over placebo years b ∈ 1930–1990 with |b − 1960| > 5, and β2(1960) > 0.
- **Kill:** F(1960) ≤ that placebo maximum.
- **Simulation:** each real child keeps its year and its number of dated grex parents. Parent ages are drawn from a gamma distribution (shape 2) with mean μ(t) = a + b₁(t − 1900) + c(t − 1900)² + δ·(t − 1960)₊, over the grid a ∈ {10, 20}, b₁ ∈ {0, 0.15}, c ∈ {0, 0.001, −0.001}, δ ∈ {0, 0.1, 0.2}, with 200 replicates per cell.
- **Open the real B only if** power (pass rate at δ = 0.2) ≥ 0.8 in every (a, b₁, c) cell, **and** the false-positive rate (pass rate at δ = 0) ≤ 0.10 in every cell, including the curved ones. Otherwise B is reported as underpowered or confounded by smooth curvature, and the real statistic stays closed.
- Clone names are absent from the scrape, so the historian's clone-string readout is unavailable (noted 2026-10-25).

## D: F1 fertility ladder by rank class
- **Outcome:** for each registered primary (two species parents), whether it was later used as a parent of any registered grex (binary). The F2/backcross rung is reported descriptively only.
- **Model:** logistic, used ~ different_rank + log1p(popularity of seed parent) + log1p(popularity of pollen parent) + decade of registration (factor). Popularity = species `primary_T1`.
- **Simulation:** real primaries keep their pair, rank class, year and parent popularities. The outcome is drawn with logit p = α + γ·different_rank + κ·(log1p pop₁ + log1p pop₂), κ ∈ {0, 0.5}, γ ∈ {0, log 0.7, log 0.5}, α set for a base rate of 0.4, 500 replicates per cell.
- **Open the real D only if** power at γ = log 0.7 (Wald p < 0.05, correct sign) ≥ 0.8 at both κ, **and** FPR at γ = 0 ≤ 0.07 at both κ. Otherwise D is reported as underpowered and the real fit stays closed.
- **Interpretation limit:** "used as a parent" mixes fertility with breeder interest. Popularity covariates absorb only the species-level part.

## E full model
Not built: the E gate returned NO-GO on readout 1 (2026-10-25, `~/orchid-data/cattleya_deeper/E_gate.json`).
