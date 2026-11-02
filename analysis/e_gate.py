"""Natural (WCVP nothospecies) vs registered primary hybrids: the E go/no-go gate.

What: maps WCVP accepted Cattleya nothospecies and register primaries (both parents species) onto
WCVP accepted species pairs, then computes the three frozen gate readouts and the 2x2 tables.
When: deciding whether the full E pair model is worth building (prereg/E_gate.md).
Output: GO needs all three readouts to pass; a NO-GO names the readout that failed.
Test: tests/test_e_gate.py (toy WCVP + toy register with planted synonyms, remakes, allopatry).
Usage: python analysis/e_gate.py HYBRIDS.csv SPECIES.csv WCVP_NAMES WCVP_DIST RANKS.csv OUT_JSON
"""

import collections
import csv
import itertools
import json
import os
import signal
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402

ABBR = {"C.": "Cattleya", "L.": "Laelia", "S.": "Sophronitis"}
csv.field_size_limit(10**8)


class Names:
    """WCVP name → accepted plant_name_id (None if absent or ambiguous)."""

    def __init__(self, rows):
        self.by_name = collections.defaultdict(set)
        self.row = {}
        for d in rows:
            self.row[d["plant_name_id"]] = d
        for d in rows:
            acc = (
                d["plant_name_id"]
                if d["taxon_status"] == "Accepted"
                else d["accepted_plant_name_id"]
            )
            if acc and acc in self.row:
                self.by_name[d["taxon_name"]].add(acc)

    def accepted(self, name):
        name = " ".join(name.split()[:2])
        hit = self.by_name.get(name, set())
        accepted_self = {a for a in hit if self.row[a]["taxon_name"] == name}
        if len(hit) == 1:
            return next(iter(hit))
        return next(iter(accepted_self)) if len(accepted_self) == 1 else None

    def name(self, pid):
        return self.row[pid]["taxon_name"]


def parse_formula(f):
    parts = [p.split() for p in f.split(" × ")]
    if len(parts) != 2 or any(len(p) < 2 or p[1] == "?" for p in parts):
        return None
    return [f"{ABBR.get(p[0], p[0])} {p[1]}" for p in parts]


def nothospecies(rows, names, status="Accepted"):
    """[(epithet, frozenset(accepted ids))], dropped count."""
    out, dropped = [], 0
    for d in rows:
        if not (
            d["genus"] == "Cattleya"
            and d["species_hybrid"] == "×"
            and d["taxon_rank"] == "Species"
            and d["taxon_status"] == status
        ):
            continue
        ps = parse_formula(d["hybrid_formula"] or "")
        ids = [names.accepted(p) for p in ps] if ps else None
        if (
            not ids
            or None in ids
            or len(set(ids)) != 2
            or any(names.name(i).split()[0] != "Cattleya" for i in ids)
        ):
            dropped += 1
            continue
        out.append((d["species"].lower(), frozenset(ids)))
    return out, dropped


def primaries(records, species, names):
    """[(epithet, grex_style, frozenset(accepted ids))] for records with two species parents."""
    idx = rq.build_index(records)
    out, dropped = [], 0
    for r in records:
        if r["synonym_of"] is not None or len(r["parents"]) != 2:
            continue
        res = [rq.resolve(p, idx, records, species) for p in r["parents"]]
        if not all(k in ("species", "species_unlisted") for k, _ in res):
            continue
        ids = [names.accepted(t) for _, t in res]
        if None in ids or len(set(ids)) != 2:
            dropped += 1
            continue
        ep = (
            rq.TAG.sub("", r["name"]).split(" ", 1)[1]
            if " " in r["name"]
            else r["name"]
        )
        out.append((ep.lower(), ep[:1].isupper(), frozenset(ids)))
    return out, dropped


def native_ranges(dist_rows):
    rng = collections.defaultdict(set)
    for d in dist_rows:
        if (
            d["introduced"] == "0"
            and d["extinct"] == "0"
            and d["location_doubtful"] == "0"
            and d["area_code_l3"]
        ):
            rng[d["plant_name_id"]].add(d["area_code_l3"])
    return rng


def gate(notho, prim, universe, ranges, rank, artificial=()):
    universe = {u for u in universe if ranges.get(u)}
    notho_ep = collections.defaultdict(set)
    for ep, pair in notho:
        notho_ep[pair].add(ep)
    crossable = {
        pair for ep, grex, pair in prim if grex and ep not in notho_ep.get(pair, ())
    }
    nat = {pair for _, pair in notho}
    nat_u = [p for p in nat if p <= universe]
    symp = lambda p: bool(set.intersection(*(ranges[i] for i in p)))
    r1 = sum(p in crossable for p in nat_u) / len(nat_u) if nat_u else float("nan")
    allo = sum(not symp(p) for p in nat_u) / len(nat_u) if nat_u else float("nan")
    cross_u = [p for p in crossable if p <= universe]
    symp_c = sum(symp(p) for p in cross_u) / len(cross_u) if cross_u else float("nan")
    art_u = [p for _, p in artificial if p <= universe]
    out = {
        "n_universe_species": len(universe),
        "n_nothospecies_in_universe": len(nat_u),
        "n_crossable_pairs": len(cross_u),
        "readout1_notho_with_independent_primary": r1,
        "readout1_pass": r1 >= 0.75,
        "readout2_notho_allopatric": allo,
        "readout2_pass": allo <= 0.20,
        "readout2_flag": 0.15 < allo <= 0.20,
        "readout3_crossable_sympatric": symp_c,
        "readout3_pass": symp_c <= 0.80,
        "artificial_allopatric": (sum(not symp(p) for p in art_u) / len(art_u))
        if art_u
        else None,
        "n_artificial_in_universe": len(art_u),
    }
    out["GO"] = out["readout1_pass"] and out["readout2_pass"] and out["readout3_pass"]
    # what the register holds for each natural-hybrid pair (descriptive; readout 1 counts "independent")
    held = collections.defaultdict(set)
    for ep, grex, pair in prim:
        held[pair].add(
            "remake"
            if grex and ep in notho_ep.get(pair, ())
            else "grex"
            if grex
            else "species_style"
        )
    kinds = collections.Counter(
        "independent"
        if p in crossable
        else "remake_under_notho_epithet"
        if "remake" in held[p]
        else "species_style_only"
        if held[p]
        else "none"
        for p in nat_u
    )
    out["readout1_breakdown"] = {
        k: kinds.get(k, 0)
        for k in (
            "independent",
            "remake_under_notho_epithet",
            "species_style_only",
            "none",
        )
    }
    tables = {}
    for a, b in itertools.combinations(sorted(universe), 2):
        p = frozenset((a, b))
        strata = ["pooled"]
        if rank.get(a) and rank.get(b):
            strata.append("same_rank" if rank[a] == rank[b] else "different_rank")
        for s in strata:
            t = tables.setdefault(s, [[0, 0], [0, 0]])
            t[0 if p in nat else 1][0 if p in crossable else 1] += 1
    out["tables"] = (
        tables  # rows: nothospecies yes/no; cols: independent primary yes/no
    )
    return out


def fisher(tables):
    from scipy.stats import fisher_exact

    res = {}
    for s, t in tables.items():
        o, p = fisher_exact(t)
        res[s] = {"table": t, "odds_ratio": o, "p": p}
    return res


def main(hybrids_csv, species_csv, names_csv, dist_csv, ranks_csv, out_json):
    rows = [
        d
        for d in csv.DictReader(open(names_csv, encoding="utf-8"), delimiter="|")
        if d["family"] == "Orchidaceae"
    ]
    names = Names(rows)
    notho, nd = nothospecies(rows, names)
    art, ad = nothospecies(rows, names, status="Artificial Hybrid")
    records, species = rq.load(hybrids_csv, species_csv)
    prim, pd = primaries(records, species, names)
    ids = {r["plant_name_id"] for r in rows}
    ranges = native_ranges(
        d
        for d in csv.DictReader(open(dist_csv, encoding="utf-8"), delimiter="|")
        if d["plant_name_id"] in ids
    )
    sp_ids = {s: names.accepted(s) for s in species}
    universe = {
        i for i in sp_ids.values() if i and names.name(i).startswith("Cattleya ")
    }
    rank = {}
    for r in csv.DictReader(open(ranks_csv, encoding="utf-8")):
        i = names.accepted(r["Binomial"])
        if i and r.get("LastSubgenus"):
            rank[i] = r["LastSubgenus"]
    res = gate(notho, prim, universe, ranges, rank, art)
    res["fisher"] = fisher(res.pop("tables"))
    res["dropped"] = {
        "nothospecies_formula": nd,
        "artificial_formula": ad,
        "primaries_unmapped": pd,
        "register_species_unmapped": sum(v is None for v in sp_ids.values()),
        "universe_without_range": len(universe) - res["n_universe_species"],
    }
    res["n_nothospecies_accepted_mapped"] = len(notho)
    res["n_primaries_mapped"] = len(prim)
    res["n_ranked_species"] = len(rank)
    with open(out_json, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "fisher"}, indent=1))
    print(json.dumps(res["fisher"]))


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(1800)
    main(*sys.argv[1:7])
