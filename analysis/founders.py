"""Founder-species composition and pedigree diversity of the register (pre-registration A).

What: expected founder-species composition of every grex (½ seed + ½ pollen, species = founders,
unresolvable parents = "unknown" mass), per-cohort founder shares, effective number of founder
species f_e = 1/Σp², and effective number of ancestors f_a (Boichard 1997 marginal contributions).
When: pedigree-diversity questions on the register (A, and generation interval for B).
Output: a composition is an expectation, not a transmitted genome: a parent grex is one clone.
Test: tests/test_founders.py (exact toy pedigrees + simulated register with planted structure).
Usage: python analysis/founders.py HYBRIDS.csv SPECIES.csv CRAWL.csv OUT_DIR
"""

import collections
import csv
import json
import os
import random
import signal
import sys

sys.path.insert(0, os.path.dirname(__file__))
import register_qc as rq  # noqa: E402

UNK = "unknown"


def binomial(ref):
    return " ".join(ref.split()[:2])


def build(records, species):
    """Per record: (year, [parent, parent]); parent = ('sp', binomial) | ('g', idx) | ('u', None)."""
    idx = rq.build_index(records)
    ped = {}
    for i, r in enumerate(records):
        if r["synonym_of"] is not None:
            continue
        ps = []
        for p in r["parents"]:
            k, t = rq.resolve(p, idx, records, species)
            for _ in range(5):  # follow synonym records to their accepted grex
                if k == "grex" and records[t]["synonym_of"] is not None:
                    k, t = rq.resolve(records[t]["synonym_of"], idx, records, species)
            if k in ("species", "species_unlisted"):
                ps.append(("sp", binomial(t)))
            elif k == "grex" and records[t]["synonym_of"] is None and t != i:
                ps.append(("g", t))
            else:
                ps.append(("u", None))
        if len(ps) == 2:
            ped[i] = (r["year"], ps)
    return ped


def topo(ped):
    kids = collections.defaultdict(list)
    indeg = {i: 0 for i in ped}
    for i, (_, ps) in ped.items():
        for k, t in ps:
            if k == "g" and t in ped:
                kids[t].append(i)
                indeg[i] += 1
    order, stack = [], [i for i, d in indeg.items() if d == 0]
    while stack:
        v = stack.pop()
        order.append(v)
        for c in kids[v]:
            indeg[c] -= 1
            if indeg[c] == 0:
                stack.append(c)
    assert len(order) == len(ped), "pedigree has a cycle"
    return order


def compositions(ped, order=None):
    comp = {}
    for i in order or topo(ped):
        out = collections.Counter()
        for k, t in ped[i][1]:
            if k == "sp":
                out[t] += 0.5
            elif k == "g" and t in comp:
                for s, v in comp[t].items():
                    out[s] += 0.5 * v
            else:
                out[UNK] += 0.5
        comp[i] = out
    return comp


def shares(comp, cohort):
    tot = collections.Counter()
    for i in cohort:
        tot.update(comp[i])
    n = len(cohort)
    return {s: v / n for s, v in tot.items()}


def f_e(sh):
    known = {s: v for s, v in sh.items() if s != UNK}
    z = sum(known.values())
    return 1.0 / sum((v / z) ** 2 for v in known.values()) if z else float("nan")


def f_a(ped, cohort, order=None, thr=0.99, max_iter=2000):
    """Effective number of ancestors by iterative marginal contributions (Boichard 1997).

    Candidates are grex ancestors of the cohort and founder species ('sp', name). Each round, cohort
    weight flows up (½ per parent) and stops at already-selected ancestors; a candidate's marginal
    is the weight reaching it (its own cohort weight excluded) × (1 − its genome fraction from
    selected ancestors). Stops at thr × the mass not lost to unresolved parents. Returns (f_a, [(ancestor, marginal), ...]);
    f_a = (Σq)² / Σq², renormalised to the explained mass as f_e is to the known mass."""
    order = order or topo(ped)
    cohort = set(cohort)
    n = len(cohort)
    anc = set(cohort)
    for i in reversed(order):
        if i in anc:
            anc.update(t for k, t in ped[i][1] if k == "g")
    sub = [i for i in order if i in anc]
    # explainable mass: everything except flow into unresolved parents ("u"), which no ancestor holds
    w0 = collections.Counter({i: 1.0 / n for i in cohort})
    lost = 0.0
    for i in reversed(sub):
        for k, t in ped[i][1]:
            if k == "g":
                w0[t] += 0.5 * w0[i]
            elif k == "u":
                lost += 0.5 * w0[i]
    known = 1.0 - lost
    sel, q = set(), []
    while len(q) < max_iter and sum(v for _, v in q) < thr * known:
        w = collections.Counter({i: 1.0 / n for i in cohort})
        for i in reversed(sub):
            if i in sel or not w[i]:
                continue
            for k, t in ped[i][1]:
                if k == "sp":
                    w[("sp", t)] += 0.5 * w[i]
                elif k == "g":
                    w[t] += 0.5 * w[i]
        g = {}
        for i in sub:
            g[i] = (
                1.0
                if i in sel
                else sum(
                    0.5
                    * (
                        float(("sp", t) in sel)
                        if k == "sp"
                        else (1.0 if t in sel else g.get(t, 0.0))
                        if k == "g"
                        else 0.0
                    )
                    for k, t in ped[i][1]
                )
            )  # a parent grex lacking parentage is not in g
        best, bv = None, 1e-12
        for a, v in w.items():
            if a in sel:
                continue
            own = 1.0 / n if a in cohort else 0.0
            m = (v - own) * (1.0 - (0.0 if isinstance(a, tuple) else g.get(a, 0.0)))
            if m > bv:
                best, bv = a, m
        if best is None:
            break
        sel.add(best)
        q.append((best, bv))
    tot = sum(v for _, v in q)
    return (tot * tot / sum(v * v for _, v in q) if q else float("nan")), q


def simulate(rng, n=2000, years=100, n_species=40, popular=True, group_shift=True):
    """Register-like pedigree with planted structure: returns (ped, truth)."""
    species = [f"S{k}" for k in range(n_species)]
    GW = 0.3 * (n_species - 2) / (0.7 * 2)
    group = species[:2]  # 2/40 = 5% of species picks; GW lifts it to 30%
    ped, uses, born = {}, collections.Counter(), []
    pop = None
    for i in range(n):
        y = int(i * years / n)

        def pick_species():
            w = [
                GW if (group_shift and y >= 60 and s in group) else 1.0 for s in species
            ]
            return ("sp", rng.choices(species, w)[0])

        def pick_grex():
            if pop is not None and rng.random() < 0.2:
                return ("g", pop)
            older = [j for j in born if ped[j][0] < y]
            if not older:
                return pick_species()
            return ("g", rng.choices(older, [uses[j] + 1 for j in older])[0])

        ps = [pick_species() if rng.random() < 0.3 else pick_grex() for _ in range(2)]
        ped[i] = (y, ps)
        for k, t in ps:
            if k == "g":
                uses[t] += 1
        born.append(i)
        if popular and pop is None and y >= 20:
            pop = i
    return ped, {"popular": pop, "group": set(group)}


def raw_contrib(ped, cohort, order):
    """Unselected expected contribution of every ancestor to the cohort (own cohort weight excluded)."""
    cohort = set(cohort)
    w = collections.Counter({i: 1.0 / len(cohort) for i in cohort})
    for i in reversed(order):
        if w[i]:
            for k, t in ped[i][1]:
                if k in ("g", "sp"):
                    w[t if k == "g" else ("sp", t)] += 0.5 * w[i]
    for i in cohort:
        w[i] -= 1.0 / len(cohort)
    return w


def edge_swap(ped, years, rng, per_edge=10):
    """Degree-preserving parent swaps between children of the same decade (pre-registration A null).

    An edge (child, slot) is swappable when the child's year is known and the parent is a species or
    a grex of known year. Two edges swap parents only if both parents are strictly older than both
    children (species count as older than everything). Each parent keeps its number of uses and
    each child keeps two parents. Returns (new_ped, accepted_swaps)."""
    new = {i: (y, list(ps)) for i, (y, ps) in ped.items()}

    def pyear(p):
        k, t = p
        if k == "sp":
            return float("-inf")
        return years.get(t) if k == "g" else None

    by_dec = collections.defaultdict(list)
    for i, (y, ps) in new.items():
        if y is None:
            continue
        for j, p in enumerate(ps):
            if pyear(p) is not None:
                by_dec[y // 10].append((i, j))
    edges = [e for es in by_dec.values() for e in es]
    ok = 0
    for _ in range(per_edge * len(edges)):
        a, ja = edges[rng.randrange(len(edges))]
        pool = by_dec[new[a][0] // 10]
        b, jb = pool[rng.randrange(len(pool))]
        if a == b:
            continue
        pa, pb = new[a][1][ja], new[b][1][jb]
        lo = min(new[a][0], new[b][0])
        if pyear(pa) < lo and pyear(pb) < lo:
            new[a][1][ja], new[b][1][jb] = pb, pa
            ok += 1
    return new, ok


def decade_stats(ped, order, dec):
    comp = compositions(ped, order)
    out = {}
    for d, cohort in dec.items():
        if len(cohort) < 20:
            continue
        fe = f_e(shares(comp, cohort))
        fa, _ = f_a(ped, cohort, order)
        rc = raw_contrib(ped, cohort, order)
        top_g = max((v for a, v in rc.items() if not isinstance(a, tuple)), default=0.0)
        out[d] = {"f_e": fe, "f_a": fa, "ratio": fe / fa, "top_grex_contrib": top_g}
    return out


def _null_draw(args):
    ped, years, dec, seed = args
    rng = random.Random(seed)
    for _ in range(5):  # a swap can close a cycle through a register date-order error: redraw
        new, ok = edge_swap(ped, years, rng)
        try:
            order = topo(new)
        except AssertionError:
            continue
        return decade_stats(new, order, dec)
    raise RuntimeError(f"null draw {seed}: 5 cyclic redraws")


def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return num / (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5


def main(hybrids_csv, species_csv, crawl_csv, out_dir):
    records, species = rq.load(hybrids_csv, species_csv)
    ped = build(records, species)
    order = topo(ped)
    comp = compositions(ped, order)
    # gate 1: second route (OrchidRoots ancestry_json)
    key = {(records[i]["name"], records[i]["year"]): i for i in ped}
    key.update(
        {(rq.TAG.sub("", records[i]["name"]), records[i]["year"]): i for i in ped}
    )
    l1 = []
    for r in csv.DictReader(open(crawl_csv, encoding="utf-8", errors="replace")):
        if not r["grex_name"].startswith("Cattleya ") or not r["ancestry_json"]:
            continue
        y = int(r["reg_year"]) if r["reg_year"].isdigit() else None
        i = key.get((rq.norm(r["grex_name"]), y))
        if i is None:
            continue
        theirs = {
            ("Cattleya " + s[3:]) if s.startswith("C. ") else s: v
            for s, v in json.loads(r["ancestry_json"]).items()
        }
        ours = {s: v for s, v in comp[i].items() if s != UNK}
        if comp[i].get(UNK, 0) > 0:
            continue
        l1.append(
            sum(abs(ours.get(s, 0) - theirs.get(s, 0)) for s in set(ours) | set(theirs))
        )
    l1.sort()
    gate1 = {
        "matched": len(l1),
        "median_L1": l1[len(l1) // 2] if l1 else None,
        "share_L1_lt_0.05": sum(x < 0.05 for x in l1) / len(l1) if l1 else None,
    }
    gate1["pass"] = bool(l1) and gate1["median_L1"] < 0.10
    # cohorts by decade
    dec = collections.defaultdict(list)
    for i in ped:
        if ped[i][0] is not None:
            dec[ped[i][0] // 10 * 10].append(i)
    dec = dict(dec)
    years = {i: r["year"] for i, r in enumerate(records) if r["synonym_of"] is None and r["year"] is not None}
    obs = decade_stats(ped, order, dec)
    rows, long, seen = [], [], set()
    for d in sorted(dec):
        sh = shares(comp, dec[d])
        for i in dec[d]:
            seen.update(s for s in comp[i] if s != UNK)
        group = lambda names: sum(sh.get(f"Cattleya {s}", 0) for s in names)
        rows.append({"decade": d, "n": len(dec[d]), "unknown": sh.get(UNK, 0.0), "f_e": f_e(sh),
                     "f_a": obs.get(d, {}).get("f_a", float("nan")),
                     "top_grex_contrib": obs.get(d, {}).get("top_grex_contrib", float("nan")),
                     "cum_species": len(seen),
                     "mini": group(("coccinea", "walkeriana", "aclandiae", "intermedia")),
                     "labiata": group(("labiata", "mossiae", "trianae", "warscewiczii", "dowiana"))})
        long += [{"decade": d, "species": s, "share": v} for s, v in sorted(sh.items(), key=lambda x: -x[1])]
    os.makedirs(out_dir, exist_ok=True)
    for name, tab in (("founders_by_decade.csv", rows), ("founder_shares_by_decade.csv", long)):
        with open(os.path.join(out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(tab[0]))
            w.writeheader()
            w.writerows(tab)
    r1 = next(r for r in rows if r["decade"] == 1950)
    r2 = next(r for r in rows if r["decade"] == 1980)
    gate3 = {"mini_1950s": r1["mini"], "mini_1980s": r2["mini"], "labiata_1950s": r1["labiata"],
             "labiata_1980s": r2["labiata"], "pass": r2["mini"] > r1["mini"] and r2["labiata"] < r1["labiata"]}
    kill = {"rho_all_decades": spearman([r["f_e"] for r in rows], [r["cum_species"] for r in rows]),
            "rho_n_ge_20": spearman(*zip(*[(r["f_e"], r["cum_species"]) for r in rows if r["n"] >= 20]))}
    kill["killed"] = kill["rho_all_decades"] > 0.9
    res = {"gate1": gate1, "gate3": gate3, "kill": kill, "n_ped": len(ped)}
    n_null = int(os.environ.get("A_NULL_DRAWS", "200"))
    if n_null:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(int(os.environ.get("A_WORKERS", "4"))) as pool:
            draws = pool.map(_null_draw, [(ped, years, dec, 1000 + k) for k in range(n_null)])
        band = {}
        for d in obs:
            for m in ("ratio", "f_e", "top_grex_contrib"):
                v = sorted(x[d][m] for x in draws)
                lo, hi = v[int(0.025 * (len(v) - 1))], v[int(0.975 * (len(v) - 1))]
                band.setdefault(str(d), {})[m] = {"obs": obs[d][m], "lo": lo, "hi": hi,
                                                  "outside": not lo <= obs[d][m] <= hi}
        res["null"] = {"draws": n_null, "band": band}
    with open(os.path.join(out_dir, "A_result.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: res[k] for k in ("gate1", "gate3", "kill", "n_ped")}))


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(3600)
    main(*sys.argv[1:5])
