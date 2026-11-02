# Pre-registration G: register QC (frozen 2026-10-25, before real counts)

Input: the 2025-03-12 orchidroots Cattleya scrape (raw, not distributed), read by `analysis/register_qc.py`.

## Resolution rule
A parent reference resolves to a grex record when its whitespace-normalised string, including any "(YYYY)" disambiguation tag, equals a grex name exactly. Synonym records are excluded when a non-synonym record has the same name. If more than one non-synonym record shares the name, the reference is ambiguous. A reference whose second word starts lowercase is a species, and resolves against the species list. Anything else is unresolved.

Prior look disclosed: a tag-stripping resolver gave 128–141 parent-after-child cases. That count is known to be mostly an artefact of the stripping. Exact matching resolved 29,623 grex and 7,553 species references, leaving 18 unresolved.

## Checks (each a flag class)
1. `date_order`: a parent grex's registration year is later than its child's.
2. `cycle`: a record lies on a directed cycle of grex→grex parentage.
3. `self_parent`: a record lists itself as a parent.
4. `ambiguous_ref`: the reference matches more than one non-synonym record.
5. `unresolved_ref`: the reference matches no grex and no species.
6. `missing_year`, `missing_parentage`, `synonym_chain` (a synonym record pointing to another synonym or to nothing): descriptive classes.

## Controls (run through the same entry point as the real data)
- Planted copy, seed 0:
  - 50 children get their year set to 1 year before a randomly chosen resolved grex parent → `date_order`.
  - 20 cycles, made by making a grex's resolved parent's first parent point back at the grex → `cycle`.
  - 100 parent references with one inner letter substituted → `unresolved_ref`.
  - 20 tags removed from references to tagged homonyms that have an untagged twin or ≥2 records → `ambiguous_ref` or wrong resolution (counted).
- Clean twin: the unplanted copy must give flags identical to the real run, and the real run is deterministic.
- Pass: recall ≥ 95% for each planted class (date_order, cycle, unresolved_ref), and zero new flags outside the planted records beyond those mechanically implied, such as descendants. The planted-tag-removal class is reported, not gated.
- Not detectable by design: a swapped or wrong parent that causes no date or cycle violation (e.g. the Pendentive → C. × dolosa class). This is stated as a limitation.

## Outputs
- `~/orchid-data/cattleya_deeper/qc_flags.tsv`: names grexes, stays outside the repo.
- `data/qc_summary.csv`: counts per class, no names.

## External judge (kill bar A4)
The OrchidRoots organizer and volunteers spot-check 50 flags drawn at random, stratified by class, blind to class. G is dead if their precision is < 70%. Until that check runs, the report is labelled "unverified by curators". Sending it needs the user's OK.

## Amendment 1 (2026-10-25, after the first real run, before any curator check)
In the first run, 108 of the 134 `unresolved_ref` flags were species-form names missing from our 141-name species list (e.g. "Cattleya dowiana var. aurea" 70, "Cattleya brevipedunculata" 34). Those are gaps in our reference list, not register errors. They now go to a descriptive class, `species_unlisted`, that is not sent to curators. The first-run flags are kept (`qc_flags_prereg_v1.tsv`), and curator precision will be reported under both the frozen and the amended definitions. Planted controls are unchanged.

## Post-hoc note (2026-10-25, after the real run, before any curator check; not pre-registered)
`analysis/register_qc_resolve.py` gives each flag a subclass from the list itself and WCVP v15. The result: 49 of 69 date flags have a benign reading (13 natural-hybrid entries whose Year is a publication date, 12 older crosses sunk as synonyms into a later-registered parent, 8 parents registered late by their own breeder, 16 gaps of 1–2 years). 39 of 196 ambiguous names have only one candidate older than the child. 11 of 18 unresolved and 38 of 52 synonym flags point to other genera. The frozen flag list, the spot-check sheet and the 70% bar are unchanged. The curator verdicts will also be reported by subclass, as a test of these readings.
