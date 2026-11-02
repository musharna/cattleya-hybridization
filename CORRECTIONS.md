# Corrections

Changes to [writeup.md](writeup.md) after it was first posted here.

## Phylogenetic signal

The species-level result in Table 1, that hybridization clusters on the phylogeny, does not hold up
and is withdrawn. The numbers in Table 1 are reproducible (lambda 0.337 regenerates exactly with the
current code), but they measure the wrong trait with a test that cannot work on this tree.

- **Wrong trait.** Each species was scored by its total registered descendants, intergeneric ones
  included. Across 100 species these sum to 443,474, many times the 19,577 *Cattleya* hybrids. A
  descendant total mostly measures how long and how heavily a species' lineage has been bred, not how
  many hybrids the species itself formed.
- **No control for breeder effort.** Species described early have had more time to be crossed.
- **The test cannot discriminate here.** On this tree, Pagel's lambda was significant in 200 of 200
  simulated data sets that had no phylogenetic signal. A placebo trait, years since the species was
  described, gives a higher lambda (0.76) than the hybrid counts do (0.46).

The re-test was written down before it was run (`analysis/prereg/C_signal_corrected.md`). The trait
is the number of primary hybrids per species, modelled as Poisson counts with years since the
species was first described as exposure and a phylogenetic random effect (brms). On the original
tree, the phylogenetic share of variation (H2) is 0.13, with a 95% credible interval of 0.0003–0.75,
so the data fit no signal and strong signal about equally well. Under the pre-registered rule, the
claim is not supported.

Because that tree is mostly unresolved, the same test was re-run on a molecular tree built from
GenBank ITS, matK and trnL-F sequences for 91 of the 100 species (pre-registered as Amendments 2
and 3). The result does not change: H2 is 0.12 (95% interval 0.0005–0.52). It agrees on 20 bootstrap
trees pooled (0.14, 0.0006–0.56), on an ITS-only tree (0.10, 0.0002–0.48) and with all 100 species
(0.12, 0.0004–0.52). Seven of those fits were refit with longer chains to reach convergence
(Amendment 6); the values barely moved. A plastid-only
tree, with sequences for only 51 species, gives the same verdict but cannot narrow the interval
(0.39, 0.002–0.99). On the combined tree the test passed none of 10 simulated data sets with no
signal and 8 of 10 with strong signal (lambda 0.6), so strong phylogenetic clustering is unlikely,
not merely undetected. At weak signal (lambda 0.3, Amendment 4) it detected it in 5 of 10
simulations, so weak clustering can neither be shown nor ruled out with these data.

Species transferred into *Cattleya* from other genera (*Laelia*, *Sophronitis* and others) have about
as many primary hybrids per year since description as the rest: a rate ratio of 0.69 (95% interval
0.47–1.03) on the molecular tree and 0.66 (0.44–1.07) on the original tree. Dropping that term gives
an H2 of 0.25 (0.02–0.60) on the molecular tree, still short of support.

An earlier version of this section counted exposure from the year of each species' current name.
For the 52 species moved into *Cattleya* in 2006–2014 that is 2006–2014 rather than their first
description (1818–1891 for most), so their exposure was about 17 years instead of more than a
century. That error produced the rate ratio of 2.5 reported for transferred species, and it shifted
every H2 above. Exposure now uses the first publication year in WCVP (Amendment 5). The verdict is
unchanged.

Affected text: the Abstract's final sentence, the Phylogenetic Signal results and discussion,
Table 1, and the Figure 3 caption.
