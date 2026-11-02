import os
import collections
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import founders as fo  # noqa: E402
import register_qc as rq  # noqa: E402

SPECIES = {"Cattleya dowiana", "Cattleya mossiae", "Cattleya coccinea", "Cattleya labiata"}


def rec(name, year, parentage):
    return rq.parse_record({"Grex Name": name, "Year": str(year), "Parentage": parentage})


def sp(s):
    return ("sp", s)


def test_compositions_half_seed_half_pollen_and_unknown_kept():
    rs = [
        rec("Cattleya A", 1900, "Cattleya dowiana × Cattleya mossiae"),
        rec("Cattleya B", 1910, "Cattleya A × Cattleya coccinea"),
        rec("Cattleya C", 1920, "Cattleya B × Cattleya Nowhere"),
    ]
    comp = fo.compositions(fo.build(rs, SPECIES))
    assert dict(comp[1]) == {"Cattleya dowiana": 0.25, "Cattleya mossiae": 0.25, "Cattleya coccinea": 0.5}
    assert comp[2][fo.UNK] == 0.5 and abs(sum(comp[2].values()) - 1) < 1e-12


def test_fe_two_for_all_x_by_y():
    ped = {i: (2000, [sp("X"), sp("Y")]) for i in range(10)}
    assert abs(fo.f_e(fo.shares(fo.compositions(ped), list(ped))) - 2.0) < 1e-12


def test_fa_one_when_cohort_descends_from_one_grex():
    ped = {0: (1900, [sp("X"), sp("Y")])}
    ped.update({i: (2000, [("g", 0), ("g", 0)]) for i in range(1, 11)})
    fa, q = fo.f_a(ped, range(1, 11))
    assert q[0][0] == 0 and abs(fa - 1.0) < 1e-12


def test_fa_hand_worked_example():
    # C1 = A x B, C2 = A x D; A, B, D each from two private species: contributions A .5, B .25, D .25
    ped = {0: (1, [sp("a1"), sp("a2")]), 1: (1, [sp("b1"), sp("b2")]), 2: (1, [sp("d1"), sp("d2")]),
           3: (2, [("g", 0), ("g", 1)]), 4: (2, [("g", 0), ("g", 2)])}
    fa, q = fo.f_a(ped, [3, 4])
    assert q[0] == (0, 0.5) and abs(fa - 1 / (0.25 + 0.0625 + 0.0625)) < 1e-12


def test_fa_marginal_discounts_genome_from_selected_ancestor():
    # C = A x B, B = A x s: A contributes .75; B's marginal is .5 x (1 - .5) = .25
    ped = {0: (1, [sp("x"), sp("y")]), 1: (2, [("g", 0), sp("s")]), 2: (3, [("g", 0), ("g", 1)])}
    fa, q = fo.f_a(ped, [2])
    assert q[0] == (0, 0.75) and abs(q[1][1] - 0.25) < 1e-12
    assert abs(fa - 1 / (0.75**2 + 0.25**2)) < 1e-12


def test_simulated_popular_parent_ranks_first():
    ped, truth = fo.simulate(random.Random(3))
    late = [i for i in ped if ped[i][0] >= 80]
    _, q = fo.f_a(ped, late, max_iter=5)
    assert q[0][0] == truth["popular"]


def _group_rise(seed, shift):
    ped, truth = fo.simulate(random.Random(seed), popular=False, group_shift=shift)
    comp = fo.compositions(ped)
    share = lambda lo, hi: sum(fo.shares(comp, [i for i in ped if lo <= ped[i][0] < hi]).get(s, 0) for s in truth["group"])
    return share(90, 100) - share(30, 60)


def test_simulated_group_shift_shows_in_shares():
    # contrast with the no-shift twin (seeds 4-9: shifted rise 0.074-0.125, twin within 0.015)
    assert _group_rise(4, True) > 0.05


def test_no_shift_twin_stays_flat():
    assert abs(_group_rise(4, False)) < 0.02


def test_fa_parent_without_parentage_is_an_ancestor_and_discounts_children():
    # grex 99 has no parentage (not a pedigree node). C = 99 x B, B = 99 x s.
    # 99 gets .5 + .25 = .75; once selected, B's marginal is .5 x (1 - .5) = .25.
    ped = {1: (2, [("g", 99), sp("s")]), 2: (3, [("g", 99), ("g", 1)])}
    fa, q = fo.f_a(ped, [2])
    assert q[0] == (99, 0.75) and abs(q[1][1] - 0.25) < 1e-12


def test_fa_unresolved_parent_mass_is_not_chased():
    # half of every cohort member is an unresolved parent: f_a stops once the other half is explained
    ped = {0: (1, [sp("x"), sp("y")])}
    ped.update({i: (2, [("g", 0), ("u", None)]) for i in range(1, 6)})
    fa, q = fo.f_a(ped, range(1, 6))
    assert q == [(0, 0.5)] and fa == 1.0


def test_edge_swap_preserves_degrees_and_strict_age():
    ped, _ = fo.simulate(random.Random(5), n=600)
    years = {i: y for i, (y, _) in ped.items()}
    new, ok = fo.edge_swap(ped, years, random.Random(1))
    uses = lambda p: sorted(collections.Counter(tuple(x) for _, ps in p.values() for x in ps).items())
    assert ok > 1000 and uses(new) == uses(ped)
    assert all(len(ps) == 2 for _, ps in new.values())
    changed = 0
    for i, (y, ps) in new.items():
        for j, p in enumerate(ps):
            if p != ped[i][1][j]:
                changed += 1
                assert p[0] == "sp" or years[p[1]] < y
    assert changed > 100
    fo.topo(new)


def test_edge_swap_never_moves_an_edge_across_decades():
    # decade 0 children use only species S0; decade 1 children only S1: swaps stay within a decade
    ped = {i: (5, [("sp", "S0"), ("sp", "S0")]) for i in range(20)}
    ped.update({i: (15, [("sp", "S1"), ("sp", "S1")]) for i in range(20, 40)})
    new, ok = fo.edge_swap(ped, {}, random.Random(2))
    assert ok > 0 and new == {i: (y, list(ps)) for i, (y, ps) in ped.items()}


def test_spearman_matches_scipy():
    from scipy.stats import spearmanr
    rng = random.Random(7)
    x = [rng.random() for _ in range(30)]
    y = [v + rng.random() for v in x]
    y[3] = y[4]
    assert abs(fo.spearman(x, y) - spearmanr(x, y).correlation) < 1e-12
