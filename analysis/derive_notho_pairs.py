"""Species pairs with a WCVP natural hybrid (accepted Cattleya nothospecies), for prereg W sensitivity (a).

What: e_gate.nothospecies on WCVP v15, mapped back to accepted binomials. Output data/nothospecies_pairs.csv
(sp1, sp2, epithet), species names only.
Usage: python -I analysis/derive_notho_pairs.py WCVP_NAMES.csv OUT.csv
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import e_gate as eg  # noqa: E402


def main(names_csv, out_csv):
    rows = [
        d
        for d in csv.DictReader(open(names_csv, encoding="utf-8"), delimiter="|")
        if d["family"] == "Orchidaceae"
    ]
    names = eg.Names(rows)
    notho, dropped = eg.nothospecies(rows, names)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sp1", "sp2", "epithet"])
        for ep, ids in sorted(notho):
            a, b = sorted(names.name(i) for i in ids)
            w.writerow([a, b, ep])
    print(
        f"nothospecies kept {len(notho)}, dropped {dropped} (E gate: 104 accepted names)"
    )


if __name__ == "__main__":
    main(*sys.argv[1:3])
