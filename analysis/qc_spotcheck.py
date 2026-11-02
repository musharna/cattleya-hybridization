"""Draw the blind curator spot-check sample from qc_flags.tsv (pre-registration G, kill bar A4).

What: up to 50 flags, stratified by class (equal share, remainder to the largest classes), shuffled,
class column withheld in the sheet; the key file maps row id -> class for scoring afterwards.
Usage: python analysis/qc_spotcheck.py QC_FLAGS.tsv OUT_DIR [N] [SEED]
"""

import csv
import random
import sys

SCORED = ("date_order", "ambiguous_ref", "unresolved_ref", "synonym_chain", "missing_parentage")


def main(flags_tsv, out_dir, n=50, seed=0):
    n, rng = int(n), random.Random(int(seed))
    rows = list(csv.DictReader(open(flags_tsv), delimiter="\t"))
    by = {c: [r for r in rows if r["class"] == c] for c in SCORED}
    k = {c: min(len(by[c]), n // len(SCORED)) for c in SCORED}
    for c in sorted(SCORED, key=lambda c: -len(by[c])):
        while sum(k.values()) < n and k[c] < len(by[c]):
            k[c] += 1
    pick = [r for c in SCORED for r in rng.sample(by[c], k[c])]
    rng.shuffle(pick)
    with open(f"{out_dir}/spotcheck_sheet.tsv", "w") as f, open(f"{out_dir}/spotcheck_key.tsv", "w") as g:
        f.write("row\tgrex\tyear\tparentage\tcurator_verdict(error/not_error/unsure)\tnote\n")
        g.write("row\tclass\tdetail\n")
        for i, r in enumerate(pick, 1):
            f.write(f"{i}\t{r['grex']}\t{r['year']}\t{r['parentage']}\t\t\n")
            g.write(f"{i}\t{r['class']}\t{r['detail']}\n")
    print({c: k[c] for c in SCORED})


if __name__ == "__main__":
    main(*sys.argv[1:5])
