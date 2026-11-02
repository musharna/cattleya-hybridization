import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import e_gate as eg  # noqa: E402
import register_qc as rq  # noqa: E402


def w(pid, name, status="Accepted", acc="", sp_h="", formula="", genus="Cattleya"):
    epithet = name.split()[-1]
    return {
        "plant_name_id": pid,
        "taxon_name": name,
        "taxon_status": status,
        "accepted_plant_name_id": acc,
        "genus": genus,
        "species_hybrid": sp_h,
        "species": epithet,
        "taxon_rank": "Species",
        "hybrid_formula": formula,
        "family": "Orchidaceae",
    }


ROWS = [
    w("1", "Cattleya aa"),
    w("2", "Cattleya bb"),
    w("3", "Cattleya cc"),
    w("4", "Cattleya dd"),
    w("5", "Laelia bb", status="Synonym", acc="2"),
    # natural hybrids: one written with a synonym and reversed order, one allopatric, one 3-parent
    w("10", "Cattleya × ab", sp_h="×", formula="L. bb × C. aa"),
    w("11", "Cattleya × cd", sp_h="×", formula="C. cc var. alba × C. dd"),
    w("12", "Cattleya × abc", sp_h="×", formula="C. aa × C. bb × C. cc"),
    w(
        "13",
        "Cattleya × art",
        status="Artificial Hybrid",
        sp_h="×",
        formula="C. aa × C. dd",
    ),
]
RANGES = {"1": {"BZL"}, "2": {"BZL"}, "3": {"PER"}, "4": {"COL"}}
SPECIES = {"Cattleya aa", "Cattleya bb", "Cattleya cc", "Cattleya dd"}


def rec(name, parentage):
    return rq.parse_record({"Grex Name": name, "Year": "1950", "Parentage": parentage})


def run(records, ranges=RANGES, rows=ROWS):
    names = eg.Names(rows)
    notho, nd = eg.nothospecies(rows, names)
    art, _ = eg.nothospecies(rows, names, status="Artificial Hybrid")
    prim, _ = eg.primaries(records, SPECIES, names)
    return eg.gate(
        notho,
        prim,
        {"1", "2", "3", "4"},
        ranges,
        {"1": "X", "2": "X", "3": "Y", "4": "Y"},
        art,
    ), nd


def test_formula_synonym_reversed_order_and_var_map_to_accepted_pairs():
    names = eg.Names(ROWS)
    notho, dropped = eg.nothospecies(ROWS, names)
    assert sorted((e, sorted(p)) for e, p in notho) == [
        ("ab", ["1", "2"]),
        ("cd", ["3", "4"]),
    ]
    assert dropped == 1  # the three-parent formula


def test_remake_under_nothospecies_epithet_is_not_independent():
    remake_only = [
        rec("Cattleya Ab", "Cattleya aa × Cattleya bb"),
        rec("Cattleya ab", "Cattleya bb × Cattleya aa"),
    ]
    g, _ = run(remake_only)
    assert g["readout1_notho_with_independent_primary"] == 0.0
    assert g["readout1_breakdown"] == {
        "independent": 0,
        "remake_under_notho_epithet": 1,
        "species_style_only": 0,
        "none": 1,
    }
    g, _ = run(
        remake_only
        + [
            rec("Cattleya Sunrise", "Cattleya aa × Cattleya bb"),
            rec("Cattleya cd", "Cattleya cc × Cattleya dd"),
        ]
    )
    assert g["readout1_notho_with_independent_primary"] == 0.5
    assert g["readout1_breakdown"] == {
        "independent": 1,
        "remake_under_notho_epithet": 0,
        "species_style_only": 1,
        "none": 0,
    }


def test_readouts_and_tables_on_planted_register():
    recs = [
        rec("Cattleya Sunrise", "Cattleya aa × Cattleya bb"),
        rec("Cattleya Dusk", "Cattleya cc × Cattleya dd"),
        rec("Cattleya Noon", "Cattleya aa × Cattleya cc"),
    ]
    g, _ = run(recs)
    assert g["readout1_notho_with_independent_primary"] == 1.0 and g["readout1_pass"]
    assert (
        g["readout2_notho_allopatric"] == 0.5 and not g["readout2_pass"]
    )  # cc (PER) x dd (COL)
    assert abs(g["readout3_crossable_sympatric"] - 1 / 3) < 1e-12 and g["readout3_pass"]
    assert g["artificial_allopatric"] == 1.0  # aa (BZL) x dd (COL)
    assert g["tables"]["pooled"] == [[2, 0], [1, 3]]
    assert g["tables"]["same_rank"] == [[2, 0], [0, 0]] and g["tables"][
        "different_rank"
    ] == [[0, 0], [1, 3]]
    assert not g["GO"]


def test_species_style_record_is_not_an_independent_primary():
    # a lowercase, species-style register name (an imported natural hybrid) is not a man-made cross
    g, _ = run([rec("Cattleya hybrida", "Cattleya aa × Cattleya bb")])
    assert (
        g["readout1_notho_with_independent_primary"] == 0.0
        and g["n_crossable_pairs"] == 0
    )


def test_power_readout_fails_when_all_crossable_pairs_sympatric():
    recs = [
        rec("Cattleya Sunrise", "Cattleya aa × Cattleya bb"),
        rec("Cattleya Dusk", "Cattleya cc × Cattleya dd"),
    ]
    everywhere = {k: {"BZL"} for k in RANGES}
    g, _ = run(recs, ranges=everywhere)
    assert g["readout3_crossable_sympatric"] == 1.0 and not g["readout3_pass"]
    assert g["readout1_pass"] and g["readout2_pass"] and not g["GO"]


def test_go_when_all_three_readouts_pass():
    recs = [
        rec("Cattleya Sunrise", "Cattleya aa × Cattleya bb"),
        rec("Cattleya Dusk", "Cattleya cc × Cattleya dd"),
        rec("Cattleya Far", "Cattleya aa × Cattleya cc"),
        rec("Cattleya Farther", "Cattleya bb × Cattleya dd"),
    ]
    ranges = dict(
        RANGES, **{"4": {"PER"}}
    )  # cd becomes sympatric; aa-cc and bb-dd stay allopatric
    g, _ = run(recs, ranges=ranges)
    assert (
        g["readout1_pass"]
        and g["readout2_notho_allopatric"] == 0.0
        and g["readout3_crossable_sympatric"] == 0.5
    )
    assert g["GO"]


def test_native_ranges_ignore_introduced_extinct_doubtful():
    d = lambda pid, a, i="0", e="0", q="0": {
        "plant_name_id": pid,
        "area_code_l3": a,
        "introduced": i,
        "extinct": e,
        "location_doubtful": q,
    }
    r = eg.native_ranges(
        [
            d("1", "BZL"),
            d("1", "FLA", i="1"),
            d("1", "PER", q="1"),
            d("1", "COL", e="1"),
        ]
    )
    assert dict(r) == {"1": {"BZL"}}
