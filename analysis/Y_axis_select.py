"""Pick the held-out trajectory-axis grexes for pre-registration Y (amendment 1, item 2).

What: up to 30 grexes per decade (1880s-2020s) with >= 2 own gallery photos that are never a parent of a
photographed grex; grexes that are not offspring in the analysis (a parent lacks a photo) come first. Their
second own photo builds the trajectory axis; Y_appearance.py drops them from the pool and offspring.
Output: OUT_JSON {decade: [[grex_pid, second_thumb_url], ...]} (local; pids only, no names).
Usage: python -I analysis/Y_axis_select.py HYBRIDS.csv SPECIES.csv PHOTO_DIR OUT_JSON
"""

import csv
import json
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import founders  # noqa: E402
import register_qc as rq  # noqa: E402

PER_DECADE = 30
SEED = 20261101


def main(hybrids_csv, species_csv, photo_dir, out_json):
    records, species = rq.load(hybrids_csv, species_csv)
    with open(hybrids_csv, encoding="utf-8", errors="replace") as f:
        links = [r["Grex Link"] for r in csv.DictReader(f)]
    assert len(links) == len(records)
    pids = [(re.search(r"/orchidaceae/(\d+)/", lk) or [None, None])[1] for lk in links]
    ped = founders.build(records, species)
    man = {}
    for line in open(os.path.join(photo_dir, "manifest.jsonl")):
        r = json.loads(line)
        man[r["grex_id"]] = r
    sp_photo = {
        os.path.splitext(f)[0].replace("_", " ")
        for f in os.listdir(os.path.join(photo_dir, "species_thumbs"))
    }
    photo = {
        i
        for i in ped
        if man.get(pids[i], {}).get("status") == "ok" and records[i]["year"] is not None
    }

    def has_photo(ref):
        k, t = ref
        return (k == "g" and t in photo) or (k == "sp" and t in sp_photo)

    parents = {t for i in photo for k, t in ped[i][1] if k == "g"}
    offspring = {i for i in photo if all(has_photo(p) for p in ped[i][1])}
    rng = random.Random(SEED)
    out = {}
    for D in range(1880, 2030, 10):
        cand = [
            i
            for i in sorted(photo)
            if D <= records[i]["year"] < D + 10
            and i not in parents
            and len(man[pids[i]].get("gallery", [])) >= 2
        ]
        rng.shuffle(cand)
        cand.sort(
            key=lambda i: i in offspring
        )  # stable: non-offspring first, random within each group
        pick = cand[:PER_DECADE]
        out[str(D)] = [
            [
                pids[i],
                man[pids[i]]["gallery"][1].replace(
                    "/static/utils/images/", "/static/utils/thumbs/", 1
                ),
            ]
            for i in pick
        ]
        print(
            D,
            len(cand),
            len(pick),
            sum(i in offspring for i in pick),
            "offspring among picks",
        )
    with open(out_json, "w") as f:
        json.dump(out, f, indent=0)
    print(sum(len(v) for v in out.values()), "axis grexes")


if __name__ == "__main__":
    main(*sys.argv[1:5])
