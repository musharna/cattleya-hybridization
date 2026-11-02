"""Flower appearance over 150 years: decade-wise selection vs transmission split (pre-registration Y + amendments 1-2).

What: from photo embeddings (orchid-clip-v8, second engine DINOv2) and register parent links, per decade D:
selection S_D = mean parent midpoint of D's offspring - mean of the photographed pool before D; transmission
T_D = offspring mean - parent midpoint mean. Projected on the trajectory axis (PC1 of decade means of held-out
second photos of axis grexes) and the mini-catt axis (held-out photos 2-3 of 4 founder species minus the labiata
group). Validity gate (AUC of own-parent vs random species-pair midpoint), mechanism check (S > 0 on the mini-catt
axis in >= 2 of 1980s-2000s, amendment 2), parent-unit bootstrap CIs, within-decade photo shuffle null, readout labels.
When: --controls first (random-vector negative control x20 and planted positive control, frozen in amendment 1);
the real run refuses to start unless PHOTO_DIR/Y_controls.json says PASS.
Output: OUT_JSON with derived numbers only (no grex names); a failing gate or control is reported, not hidden.
Test: tests/test_Y_appearance.py (planted registers; each test seen to fail under mutate-run).
Usage: python -I analysis/Y_appearance.py HYBRIDS.csv SPECIES.csv PHOTO_DIR OUT_JSON [--controls]
"""

import csv
import json
import os
import re
import signal
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import founders  # noqa: E402
import register_qc as rq  # noqa: E402

DECADES = list(range(1880, 2030, 10))
MIN_N = 20
FOUNDERS = ["coccinea", "walkeriana", "aclandiae", "intermedia"]
LABIATA = ["labiata", "mossiae", "trianae", "warscewiczii", "dowiana"]
MECH_DECADES = [1980, 1990, 2000]  # amendment 2 (was 1960s-80s)
PHOTO_DIRS = (
    ("grex", "thumbs"),
    ("species", "species_thumbs"),
    ("axis", "axis_heldout"),
    ("axisgrex", "axis_grex"),
)


def gate_auc(z_off, z_pa, z_pb, z_species, rng, reps=200):
    """Fraction of primaries whose own parents' midpoint is closer (cosine) than a random species pair's."""
    own = _cos(z_off, (z_pa + z_pb) / 2)
    aucs = []
    n = len(z_species)
    for _ in range(reps):
        a = rng.integers(0, n, len(z_off))
        b = (a + rng.integers(1, n, len(z_off))) % n  # distinct species
        rnd = _cos(z_off, (z_species[a] + z_species[b]) / 2)
        aucs.append(float(np.mean((own > rnd) + 0.5 * (own == rnd))))
    return (
        float(np.mean(aucs)),
        float(np.percentile(aucs, 2.5)),
        float(np.percentile(aucs, 97.5)),
    )


def _cos(x, y):
    return np.sum(x * y, 1) / (np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1))


def split(units, axis, rng, n_boot=1000):
    """Per-decade S, T on axis with parent-unit bootstrap CIs (amendment 1, item 3).

    units: per decade dict(pool=(n,d), off=(n,d), pa=(n,d), pb=(n,d), ka=[key], kb=[key]); ka/kb identify parent
    photos. A replicate resamples the decade's distinct parent photos with replacement and weights each offspring
    by the product of its two parents' counts.
    """
    out = {}
    for D, u in units.items():
        if len(u["off"]) == 0:
            continue
        pool = float(u["pool"].mean(0) @ axis)
        o, m = u["off"] @ axis, ((u["pa"] + u["pb"]) / 2) @ axis
        S, T = m.mean() - pool, o.mean() - m.mean()
        keys = sorted(set(u["ka"]) | set(u["kb"]), key=repr)
        pos = {k: j for j, k in enumerate(keys)}
        ia = np.array([pos[k] for k in u["ka"]])
        ib = np.array([pos[k] for k in u["kb"]])
        draws = rng.integers(0, len(keys), (n_boot, len(keys)))
        C = np.zeros((n_boot, len(keys)))
        np.add.at(C, (np.repeat(np.arange(n_boot), len(keys)), draws.ravel()), 1)
        W = C[:, ia] * C[:, ib]
        ok = W.sum(1) > 0
        W = W[ok] / W[ok].sum(1, keepdims=True)
        bm, bo = W @ m, W @ o
        bs, bt = bm - pool, bo - bm
        mv = ((u["pa"] + u["pb"]) / 2).mean(0)
        out[D] = {
            "n": len(o),
            "n_parent_photos": len(keys),
            "n_boot_used": int(ok.sum()),
            "S": float(S),
            "S_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "T": float(T),
            "T_ci": [float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))],
            "S_norm": float(np.linalg.norm(mv - u["pool"].mean(0))),
            "T_norm": float(np.linalg.norm(u["off"].mean(0) - mv)),
        }
    return out


def label(d):
    if d["n"] < MIN_N:
        return "unlabelled (n < 20)"
    s_ex = d["S_ci"][0] > 0 or d["S_ci"][1] < 0
    t_ex = d["T_ci"][0] > 0 or d["T_ci"][1] < 0
    if abs(d["S"]) >= 2 * abs(d["T"]) and s_ex:
        return "selection-led"
    if abs(d["T"]) >= 2 * abs(d["S"]) and t_ex:
        return "transmission-led"
    if s_ex and t_ex:
        return "both"
    if not s_ex and not t_ex:
        return "no clear change"
    return "unresolved"  # amendment 1, item 4: one CI excludes 0 but the 2x rule is not met


def mechanism(mc):
    hits = [
        D
        for D in MECH_DECADES
        if D in mc and mc[D]["n"] >= MIN_N and mc[D]["S_ci"][0] > 0
    ]
    return {"decades_S_above_0": hits, "PASS": len(hits) >= 2}


def trajectory_axis(dec_means):
    """PC1 of decade means, signed so the last decade projects above the first."""
    X = np.stack(list(dec_means.values()))
    Xc = X - X.mean(0)
    v = np.linalg.svd(Xc, full_matrices=False)[2][0]
    if (X[-1] - X[0]) @ v < 0:
        v = -v
    return v / np.linalg.norm(v)


def shuffle_null(units_idx, proj, decade_of, rng, reps=500):
    """Shuffle grex photos within each decade (species photos fixed); return S per decade per rep on the axis.

    units_idx: per decade dict(pool_mean=float, pa=[(kind, key)], pb=[(kind, key)]); proj: dict key -> scalar
    for grexes ('g', idx) and species ('s', name). decade_of: grex idx -> decade.
    """
    by_dec = {}
    for g, D in decade_of.items():
        by_dec.setdefault(D, []).append(g)
    out = {D: [] for D in units_idx}
    for _ in range(reps):
        perm = {}
        for D, gs in by_dec.items():
            for a, b in zip(gs, rng.permutation(gs)):
                perm[a] = b
        for D, u in units_idx.items():
            m = np.mean(
                [
                    (_p(proj, perm, pa) + _p(proj, perm, pb)) / 2
                    for pa, pb in zip(u["pa"], u["pb"])
                ]
            )
            out[D].append(m - u["pool_mean"])
    return out


def _p(proj, perm, ref):
    kind, key = ref
    return proj[("g", perm[key])] if kind == "g" else proj[("s", key)]


def photo_keys(photo_dir):
    """(kind, id) of every photo, listed the way embed_thumbs.py lists them."""
    out = []
    for kind, sub in PHOTO_DIRS:
        d = os.path.join(photo_dir, sub)
        if os.path.isdir(d):
            out += [
                (kind, os.path.splitext(f)[0])
                for f in sorted(os.listdir(d))
                if not f.endswith(".jsonl")
            ]
    return out


def load_emb(photo_dir, engine):
    z = np.load(os.path.join(photo_dir, f"emb_{engine}.npz"), allow_pickle=False)
    return {(str(k), str(i)): e for k, i, e in zip(z["kind"], z["ids"], z["emb"])}


def build_units(records, ped, pids, sp_year, emb, axis_pids):
    """Per decade pool / offspring arrays; axis grexes are dropped from pool and offspring (amendment 1, item 2)."""
    gz = {}
    for i in ped:
        e = emb.get(("grex", pids[i]))
        if (
            e is not None
            and records[i]["year"] is not None
            and pids[i] not in axis_pids
        ):
            gz[i] = e
    sz = {}
    for name in sp_year:
        e = emb.get(("species", name.replace(" ", "_")))
        if e is not None:
            sz[name] = e

    def par(ref):
        k, t = ref
        if k == "g" and t in gz:
            return ("g", t), gz[t]
        if k == "sp" and t in sz:
            return ("s", t), sz[t]
        return None, None

    units, idx = {}, {}
    for D in DECADES:
        pool = [gz[i] for i in gz if records[i]["year"] < D] + [
            sz[s] for s in sz if sp_year[s] is not None and sp_year[s] < D
        ]
        off, pa, pb, ka, kb = [], [], [], [], []
        for i in gz:
            if not D <= records[i]["year"] < D + 10:
                continue
            (ra, a), (rb, b) = par(ped[i][1][0]), par(ped[i][1][1])
            if a is None or b is None:
                continue
            off.append(gz[i])
            pa.append(a)
            pb.append(b)
            ka.append(ra)
            kb.append(rb)
        units[D] = {
            "pool": np.array(pool),
            "off": np.array(off),
            "pa": np.array(pa),
            "pb": np.array(pb),
            "ka": ka,
            "kb": kb,
        }
        idx[D] = {"pa": ka, "pb": kb}
    return units, idx, gz, sz


def axes(records, emb, pid_idx, axis_pids):
    """Trajectory axis from axis grexes' second photos; mini-catt axis from held-out species photos 2-3."""
    dec = {}
    for pid in axis_pids:
        e = emb.get(("axisgrex", pid))
        if e is not None:
            dec.setdefault(records[pid_idx[pid]]["year"] // 10 * 10, []).append(e)
    dec_means = {D: np.mean(dec[D], 0) for D in DECADES if D in dec}
    traj = trajectory_axis(dec_means)

    def held(s):
        v = [
            emb[("axis", f"Cattleya_{s}__{k}")]
            for k in (1, 2)
            if ("axis", f"Cattleya_{s}__{k}") in emb
        ]
        assert v, f"no held-out photo for Cattleya {s}"
        return np.mean(v, 0)

    mini = np.mean([held(s) for s in FOUNDERS], 0) - np.mean(
        [held(s) for s in LABIATA], 0
    )
    return traj, mini / np.linalg.norm(mini), {D: len(v) for D, v in dec.items()}


def run_engine(records, ped, pids, sp_year, emb, axis_pids, seed, extra_axes=None):
    rng = np.random.default_rng(seed)
    pid_idx = {p: i for i, p in enumerate(pids) if p is not None}
    units, idx, gz, sz = build_units(records, ped, pids, sp_year, emb, axis_pids)
    prim = [i for i in gz if all(k == "sp" and t in sz for k, t in ped[i][1])]
    zs = np.stack(list(sz.values()))
    zo = np.stack([gz[i] for i in prim])
    za = np.stack([sz[ped[i][1][0][1]] for i in prim])
    zb = np.stack([sz[ped[i][1][1][1]] for i in prim])
    auc = gate_auc(zo, za, zb, zs, rng)
    traj, mini, n_axis = axes(records, emb, pid_idx, axis_pids)
    res_t = split(units, traj, rng)
    res_m = split(units, mini, rng)
    proj = {("g", i): float(e @ traj) for i, e in gz.items()}
    proj.update({("s", s): float(e @ traj) for s, e in sz.items()})
    decade_of = {i: records[i]["year"] // 10 * 10 for i in gz}
    uidx = {
        D: {
            "pool_mean": float(units[D]["pool"].mean(0) @ traj),
            "pa": idx[D]["pa"],
            "pb": idx[D]["pb"],
        }
        for D in res_t
    }
    null = shuffle_null(uidx, proj, decade_of, rng)
    for D, d in res_t.items():
        nl = np.array(null[D])
        d["null_p"] = float(
            (np.sum(np.abs(nl - nl.mean()) >= abs(d["S"] - nl.mean())) + 1)
            / (len(nl) + 1)
        )
        d["null_mean"] = float(nl.mean())
        d["label"] = label(d)
    out = {
        "n_grex_photos_in_pedigree": len(gz),
        "n_species_photos": len(sz),
        "n_axis_grex_photos_by_decade": {str(D): n for D, n in sorted(n_axis.items())},
        "gate": {
            "n_primaries": len(prim),
            "auc": auc[0],
            "auc_95": auc[1:],
            "PASS": auc[0] >= 0.70,
        },
        "mechanism": mechanism(res_m),
        "trajectory": {str(D): d for D, d in res_t.items()},
        "minicatt": {str(D): d for D, d in res_m.items()},
        "overall": {
            "sum_S_all": float(sum(d["S"] for d in res_t.values())),
            "sum_T_all": float(sum(d["T"] for d in res_t.values())),
            "sum_S_labelled": float(
                sum(d["S"] for d in res_t.values() if d["n"] >= MIN_N)
            ),
            "sum_T_labelled": float(
                sum(d["T"] for d in res_t.values() if d["n"] >= MIN_N)
            ),
        },
        "axis_cos_traj_minicatt": float(traj @ mini),
    }
    for name, v in (extra_axes or {}).items():
        out[name] = {str(D): d for D, d in split(units, v, rng).items()}
    return out


def synthetic(keys, records, ped, pids, rng, planted, d=256):
    """Photo embeddings for the controls. planted=False: independent random unit vectors. planted=True: species
    vectors random, founders +0.5u, labiata group -0.5u; grex = parents' midpoint + noise (+0.3w for 1990s grexes);
    every photo = its item's vector + independent noise, L2-normalized. Returns (emb, {'u': u, 'w': w})."""
    unit = lambda x: x / np.linalg.norm(x, axis=-1, keepdims=True)  # noqa: E731
    if not planted:
        return {k: unit(rng.normal(size=d)) for k in keys}, {}
    u, w = np.eye(d)[0], np.eye(d)[1]
    sp_true = {}

    def sp(name):
        if name not in sp_true:
            v = rng.normal(size=d) / np.sqrt(d)
            ep = name.split()[1] if len(name.split()) > 1 else ""
            v += 0.5 * u * ((ep in FOUNDERS) - (ep in LABIATA))
            sp_true[name] = v
        return sp_true[name]

    g_true = {}
    for i in founders.topo(ped):
        ps = []
        for k, t in ped[i][1]:
            if k == "sp":
                ps.append(sp(t))
            elif k == "g" and t in g_true:
                ps.append(g_true[t])
            else:
                ps.append(rng.normal(size=d) / np.sqrt(d))
        v = (ps[0] + ps[1]) / 2 + 0.5 * rng.normal(size=d) / np.sqrt(d)
        y = records[i]["year"]
        if y is not None and 1990 <= y < 2000:
            v = v + 0.3 * w
        g_true[i] = v
    pid_idx = {p: i for i, p in enumerate(pids) if p is not None}
    emb = {}
    for kind, ident in keys:
        if kind in ("grex", "axisgrex"):
            i = pid_idx.get(ident)
            base = g_true.get(i) if i is not None else None
        else:
            base = sp(ident.split("__")[0].replace("_", " "))
        if base is None:
            base = rng.normal(size=d) / np.sqrt(d)
        emb[(kind, ident)] = unit(base + 0.5 * rng.normal(size=d) / np.sqrt(d))
    return emb, {"u": u, "w": w}


def controls(records, ped, pids, sp_year, keys, axis_pids, n_neg=20):
    """Frozen amendment-1 controls through the production path."""
    neg = []
    for s in range(n_neg):
        emb, _ = synthetic(
            keys, records, ped, pids, np.random.default_rng(1000 + s), planted=False
        )
        r = run_engine(records, ped, pids, sp_year, emb, axis_pids, seed=s)
        lab = [d["label"] for d in r["trajectory"].values() if d["n"] >= MIN_N]
        neg.append(
            {
                "gate_auc": r["gate"]["auc"],
                "gate_pass": r["gate"]["PASS"],
                "mech_pass": r["mechanism"]["PASS"],
                "labelled": len(lab),
                "false_labels": sum(x != "no clear change" for x in lab),
            }
        )
        print(f"neg {s}: {neg[-1]}", flush=True)
    emb, dirs = synthetic(
        keys, records, ped, pids, np.random.default_rng(7), planted=True
    )
    pos = run_engine(
        records,
        ped,
        pids,
        sp_year,
        emb,
        axis_pids,
        seed=99,
        extra_axes={"w": dirs["w"]},
    )
    t90 = pos["w"]["1990"]
    cells = sum(x["labelled"] for x in neg)
    false = sum(x["false_labels"] for x in neg)
    res = {
        "negative": {
            "runs": neg,
            "gate_fail_all": all(not x["gate_pass"] for x in neg),
            "mech_passes": sum(x["mech_pass"] for x in neg),
            "false_label_rate": false / cells if cells else None,
            "labelled_cells": cells,
        },
        "positive": {
            "gate_auc": pos["gate"]["auc"],
            "mechanism": pos["mechanism"],
            "T_w_1990": {"T": t90["T"], "T_ci": t90["T_ci"], "n": t90["n"]},
        },
    }
    res["negative"]["PASS"] = (
        res["negative"]["gate_fail_all"]
        and res["negative"]["mech_passes"] <= 2
        and cells > 0
        and false / cells <= 0.10
    )
    res["positive"]["PASS"] = pos["mechanism"]["PASS"] and t90["T_ci"][0] > 0
    res["PASS"] = res["negative"]["PASS"] and res["positive"]["PASS"]
    return res


def load_all(hybrids_csv, species_csv, photo_dir):
    records, species = rq.load(hybrids_csv, species_csv)
    with open(hybrids_csv, encoding="utf-8", errors="replace") as f:
        links = [r["Grex Link"] for r in csv.DictReader(f)]
    assert len(links) == len(records)
    pids = [(re.search(r"/orchidaceae/(\d+)/", lk) or [None, None])[1] for lk in links]
    with open(species_csv, encoding="utf-8", errors="replace") as f:
        sp_year = {
            rq.norm(r["Binomial"]): (
                int(r["Year"]) if r["Year"].strip().isdigit() else None
            )
            for r in csv.DictReader(f)
        }
    ped = founders.build(records, species)
    # species parents outside the species list (natural hybrids etc.) have photos but no description year:
    # they enter midpoints, not the pool
    for _, (_, ps) in ped.items():
        for k, t in ps:
            if k == "sp":
                sp_year.setdefault(t, None)
    sel = json.load(open(os.path.join(photo_dir, "axis_grexes.json")))
    axis_pids = {pid for v in sel.values() for pid, _ in v}
    return records, ped, pids, sp_year, axis_pids


def main(hybrids_csv, species_csv, photo_dir, out_json, mode=None):
    records, ped, pids, sp_year, axis_pids = load_all(
        hybrids_csv, species_csv, photo_dir
    )
    ctrl_f = os.path.join(photo_dir, "Y_controls.json")
    if mode == "--controls":
        res = controls(records, ped, pids, sp_year, photo_keys(photo_dir), axis_pids)
        for f in (ctrl_f, out_json):
            with open(f, "w") as fo:
                json.dump(res, fo, indent=1)
        print(
            f"negative PASS={res['negative']['PASS']} (gate fails all={res['negative']['gate_fail_all']}, "
            f"mech passes {res['negative']['mech_passes']}/20, false-label rate {res['negative']['false_label_rate']:.3f} "
            f"of {res['negative']['labelled_cells']}); positive PASS={res['positive']['PASS']} "
            f"(mechanism {res['positive']['mechanism']}, 1990s T_w {res['positive']['T_w_1990']}); PASS={res['PASS']}"
        )
        return
    assert os.path.exists(ctrl_f) and json.load(open(ctrl_f))["PASS"], (
        "controls have not passed: run --controls first"
    )
    out = {}
    for engine, seed in (("v8", 1), ("dino", 2)):
        out[engine] = run_engine(
            records, ped, pids, sp_year, load_emb(photo_dir, engine), axis_pids, seed
        )
    out["engine_dependent_decades"] = [
        D
        for D in out["v8"]["trajectory"]
        if D in out["dino"]["trajectory"]
        and np.sign(out["v8"]["trajectory"][D]["S"])
        != np.sign(out["dino"]["trajectory"][D]["S"])
    ]
    with open(out_json, "w") as f:
        json.dump(out, f, indent=1)
    for e in ("v8", "dino"):
        r = out[e]
        print(
            f"{e}: gate AUC {r['gate']['auc']:.3f} (n={r['gate']['n_primaries']}) PASS={r['gate']['PASS']}; mechanism {r['mechanism']}"
        )
        for D, d in r["trajectory"].items():
            print(
                f"  {D}s n={d['n']:4d} S={d['S']:+.4f} [{d['S_ci'][0]:+.4f},{d['S_ci'][1]:+.4f}] "
                f"T={d['T']:+.4f} [{d['T_ci'][0]:+.4f},{d['T_ci'][1]:+.4f}] null_p={d['null_p']:.3f} {d['label']}"
            )
    print("engine-dependent decades:", out["engine_dependent_decades"])


if __name__ == "__main__":
    signal.signal(
        signal.SIGALRM,
        lambda *_: (sys.stderr.write("aborting: walltime guard\n"), sys.exit(2)),
    )
    signal.alarm(5400)
    main(*sys.argv[1:6])
