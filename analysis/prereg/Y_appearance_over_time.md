# Pre-registration Y: how did the look of registered Cattleya change over 150 years, and was it parent choice or offspring change? (2026-10-29)

Written while the photos are being downloaded, before any image has been embedded or looked at.

**Question.** Where in appearance space did registered *Cattleya* grexes move, decade by decade? For each decade, did the move come from **selection** (which grexes breeders chose as parents) or from **transmission** (offspring looking different from their parents' midpoint)?

**Prior art** (project memory, `cattleya_expansion_prior_art_2026-10-29`).
- **Shape change over time:** dated specimens have been used for St Bernard skulls over 120 years (Drake & Klingenberg 2008) and for pig skulls (Haruda 2025).
- **Price split:** applied to 23,000 novels without parent links (Sobchuk & Beheim 2025).
- **Open part:** a register's explicit two-parent links give a real selection term. The transmission term is close to the off-chord work in the researcher's own papers, so the selection term and the 150-year trajectory are the new parts.

**Data.**
- **Grex photos:** the first own-gallery photo (300 px thumbnail) of each grex with photos in the March 2025 scrape, fetched politely from OrchidRoots (`~/orchid-data/cattleya_deeper/photos/`, not redistributed).
- **Species-parent photos:** the first gallery photo of each species parent, fetched the same way after the grex run. All photos then come from one source.
- **Embeddings:** orchid-clip-v8 image embeddings (`musharna/orchid-clip-v8`, L2-normalized). Second engine: DINOv2 ViT-B/14 CLS embeddings.
- **Links:** parent links resolved with `register_qc`. A parent counts only if it has a photo.

**Validity gate (run first).** Among photographed primaries whose two parent species both have photos, the cosine similarity between the primary and its parents' midpoint must beat the midpoint of a random photographed species pair in at least 70% of comparisons (AUC ≥ 0.70, one random pair per primary, 200 repeats). If v8 fails this, it is not used for Y. If DINOv2 also fails, Y stops.

**Quantities, per decade D (1880s–2020s).**
- **Pool:** photographed grexes and species registered or described before D, with mean z̄_pool.
- **Offspring:** photographed grexes registered in D whose two parents both have photos. ō_D is their mean, and m̄_D is the mean of their parents' midpoints.
- **Selection** S_D = m̄_D − z̄_pool. **Transmission** T_D = ō_D − m̄_D. Total Δ_D = ō_D − z̄_pool = S_D + T_D.
- **Projection axes:**
  - (1) the trajectory axis, the first principal component of the decade means of all photographed grexes;
  - (2) the mini-catt axis, the mean of the 4 founder species' photos (*coccinea*, *walkeriana*, *aclandiae*, *intermedia*) minus the mean of the labiata group (*labiata*, *mossiae*, *trianae*, *warscewiczii*, *dowiana*).
  - Norms of S_D and T_D are reported too.
- Decades with fewer than 20 qualifying offspring are shown but not labelled.
- **Uncertainty:** a bootstrap over offspring grexes (1,000 resamples) gives 95% CIs.

**Null.** Shuffle photos among grexes within each decade (500 times). This breaks parent links but keeps the time trend. Real S_D on the trajectory axis is compared with the shuffled distribution.

**Mechanism check (known answer).** Breeders took up the mini-catt founders after 1960 (analysis A: 1% → 16% of ancestry). On the mini-catt axis, S_D must be positive, with a CI above 0, in at least 2 of the decades 1960s, 1970s and 1980s. If not, the method cannot find a selection we know happened, and S/T labels are not read.

**Readout, per labelled decade, on the trajectory axis.**
- **selection-led:** |S| ≥ 2|T| and S's CI excludes 0.
- **transmission-led:** |T| ≥ 2|S| and T's CI excludes 0.
- **both:** both CIs exclude 0 and neither dominates.
- **no clear change:** both CIs include 0.

Plus one sentence for 1880s→2020s overall: the summed S versus the summed T.

**Second engine.** The same quantities from DINOv2. Labels are reported side by side; a decade whose S sign differs between engines is marked **engine-dependent**.

**Known limits, stated up front.** Photographed grexes are a selected sample: better-known and more often used grexes are more likely to have photos. Photos vary in lighting and framing. Thumbnails are 300 px. Results describe the photographed register, not every grex.

**Figures.**
- (1) Decade means of all photographed grexes in the first two principal components, with a few small thumbnails along the path.
- (2) S_D and T_D per decade on the trajectory axis and the mini-catt axis, with CIs.

`docs/v2_deeper.md` changes only after the user reviews the result.

## Amendment 1 (2026-11-01, before any real embedding exists)

**Why.** A plumbing run on random unit vectors for every photo (real pedigree, real photo set) passed the gate correctly (AUC 0.50 v8 / 0.47 DINOv2) but showed three defects in the frozen design:
- the mechanism check passed on noise (S > 0 in the 1960s, 70s and 80s), because the mini-catt axis was built from the same species photos that enter parent midpoints;
- 13 of 15 decades got a label from noise, because the offspring bootstrap treats offspring as independent while a few parents are reused hundreds of times;
- the trajectory axis was fitted to the photos it then measures (1880s T = −0.17, CI −0.19 to −0.16, from noise).

**Changes.**
1. **Mini-catt axis from held-out photos:** gallery photos 2 and 3 of each of the 9 axis species. Photo 1 stays the parent photo.
2. **Trajectory axis from held-out photos and grexes:**
   - The axis is PC1 of decade means of the *second* own-gallery photo of up to 30 grexes per decade (1880s–2020s), signed so the 2020s project above the 1880s.
   - Axis grexes are drawn (seed 20261101) from grexes with ≥ 2 own photos that are never a parent of a photographed grex. Grexes that are not offspring in the analysis come first.
   - Axis grexes are dropped from the pool and from offspring.
3. **Uncertainty: parent-unit bootstrap.**
   - The distinct parent photos of decade D are resampled with replacement (1,000 times).
   - Each offspring is weighted by the product of its two parents' resample counts.
   - This gives the S and T CIs used by the readout rule, and it replaces the offspring bootstrap.
   - The within-decade shuffle null is unchanged.
4. **Uncovered readout case:** one CI excludes 0 but the 2× rule fails. This case is labelled **unresolved**.
5. **Controls that must pass before the real embeddings are read** (same code path, real pedigree, real photo set):
   - **Negative control (random unit vectors, 20 seeds):**
     - the gate fails in every seed;
     - the mechanism check passes in ≤ 2 of 20;
     - among labelled decades, labels other than "no clear change" are ≤ 10% of decade × seed cells.
   - **Positive control (planted):**
     - Species vectors are random. Founder species are shifted +0.5 along a planted direction u, and labiata-group species −0.5.
     - Every photo is its item's vector plus independent noise.
     - Grex = parents' midpoint + noise, plus a planted shift of +0.3 along a second direction w in 1990s offspring only.
     - Pass requires the mechanism check to pass, and the 1990s T on the w-projection to have a CI above 0.
   - **Outcome:** if either control fails, Y stops and is reported as such.

## Amendment 2 (2026-11-01, before any real embedding exists)

**Why.** A dry run of the amendment-1 controls (real pedigree; the axis-grex photos were still downloading, so their IDs came from the selection file) failed one control. The negative control passed: 2 seeds, gate failed in both, mechanism check 0 of 2, false-label rate 1 of 30. The positive control failed the mechanism check. Its planted selection is +0.5u per founder species, measured on the true planted direction with no axis error:
- 1960s: S +0.030, CI −0.059 to +0.115;
- 1970s: S +0.036, CI −0.079 to +0.152;
- 1980s: S +0.113, CI +0.010 to +0.216.

This is not a code fault. It is the register's timing: founder species were 3.9% of ancestry in the 1960s and 6.5% in the 1970s, then 16% in the 1980s and 22% in the 1990s (analysis A, `data/founder_shares_by_decade.csv`). The frozen window sat where the known selection is too weak to resolve with parent-unit uncertainty.

**Change.** The mechanism-check window becomes the 1980s, 1990s and 2000s. The rule stays S > 0 with CI above 0 in at least 2 of the 3, on the held-out mini-catt axis. The window is chosen from analysis A's founder shares (≥ 16%), not from the control's output. The planted effect sizes, the negative-control thresholds and every other rule are unchanged. The amendment-1 controls are rerun in full (20 negative seeds) on the complete photo set before the real embeddings are read. If either control fails then, Y stops.
