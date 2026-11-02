import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import register_qc as rq  # noqa: E402

SPECIES = {
    "Cattleya dowiana",
    "Cattleya mossiae",
    "Cattleya coccinea",
    "Cattleya labiata",
}


def rec(name, year, parentage):
    return rq.parse_record(
        {"Grex Name": name, "Year": str(year), "Parentage": parentage}
    )


def flags_of(records, cls):
    return {i for i, c, _ in rq.qc(records, SPECIES) if c == cls}


def test_year_tag_resolves_to_tagged_homonym_not_untagged():
    rs = [
        rec("Cattleya Doris (1894)", 1894, "Cattleya dowiana × Cattleya mossiae"),
        rec("Cattleya Doris", 1950, "Cattleya coccinea × Cattleya labiata"),
        rec("Cattleya Child", 1900, "Cattleya Doris (1894) × Cattleya dowiana"),
    ]
    idx = rq.build_index(rs)
    kind, tgt = rq.resolve("Cattleya Doris (1894)", idx, rs, SPECIES)
    assert (kind, tgt) == ("grex", 0)
    # the tagged parent (1894) precedes its 1900 child: no date-order flag
    assert flags_of(rs, "date_order") == set()


def test_date_order_flag_and_clean_twin():
    rs = [
        rec("Cattleya Late", 1950, "Cattleya dowiana × Cattleya mossiae"),
        rec("Cattleya Early", 1920, "Cattleya Late × Cattleya coccinea"),
    ]
    assert flags_of(rs, "date_order") == {1}
    ok = [rs[0], rec("Cattleya Early", 1960, "Cattleya Late × Cattleya coccinea")]
    assert flags_of(ok, "date_order") == set()


def test_synonym_record_excluded_when_nonsynonym_shares_name():
    rs = [
        rec("Cattleya Wel", 1909, "Cattleya labiata × Cattleya coccinea"),
        rec("Cattleya Wel", 1907, "Synonym of Cattleya Callisto"),
        rec("Cattleya Kid", 1930, "Cattleya Wel × Cattleya dowiana"),
    ]
    idx = rq.build_index(rs)
    assert rq.resolve("Cattleya Wel", idx, rs, SPECIES) == ("grex", 0)


def test_two_untagged_nonsynonyms_are_ambiguous():
    rs = [
        rec("Cattleya Ava", 1922, "Cattleya dowiana × Cattleya mossiae"),
        rec("Cattleya Ava", 1955, "Cattleya coccinea × Cattleya labiata"),
        rec("Cattleya Kid", 1960, "Cattleya Ava × Cattleya dowiana"),
    ]
    assert flags_of(rs, "ambiguous_ref") == {2}


def test_cycle_detected_and_acyclic_twin_clean():
    rs = [
        rec("Cattleya A", 1950, "Cattleya B × Cattleya dowiana"),
        rec("Cattleya B", 1940, "Cattleya A × Cattleya mossiae"),
        rec("Cattleya C", 1960, "Cattleya A × Cattleya B"),
    ]
    assert flags_of(rs, "cycle") == {0, 1}
    ok = [
        rec("Cattleya A", 1950, "Cattleya B × Cattleya dowiana"),
        rec("Cattleya B", 1940, "Cattleya labiata × Cattleya mossiae"),
        rs[2],
    ]
    assert flags_of(ok, "cycle") == set()


def test_self_parent_and_unresolved_and_whitespace():
    rs = [
        rec("Cattleya Zoe (1934)", 1934, "Cattleya dowiana × Cattleya mossiae"),
        rec("Cattleya Me", 1950, "Cattleya Me × Cattleya dowiana"),
        rec("Cattleya X", 1960, "Cattleya \tZoe (1934) × Catyclia Whitei"),
    ]
    assert flags_of(rs, "self_parent") == {1}
    idx = rq.build_index(rs)
    assert rq.resolve("Cattleya \tZoe (1934)", idx, rs, SPECIES) == ("grex", 0)
    assert flags_of(rs, "unresolved_ref") == {2}


def test_species_reference_resolves():
    rs = [rec("Cattleya P", 1900, "Cattleya dowiana × Cattleya mossiae")]
    idx = rq.build_index(rs)
    assert rq.resolve("Cattleya dowiana", idx, rs, SPECIES) == (
        "species",
        "Cattleya dowiana",
    )
    assert rq.qc(rs, SPECIES) == []


def test_planting_recall_on_synthetic_register():
    rng = random.Random(0)
    rs = [
        rec(f"Cattleya G{i}", 1900 + i, "Cattleya dowiana × Cattleya mossiae")
        for i in range(60)
    ]
    for i in range(60, 400):
        a, b = rng.sample(range(i), 2)
        rs.append(rec(f"Cattleya G{i}", 1900 + i, f"Cattleya G{a} × Cattleya G{b}"))
    assert rq.qc(rs, SPECIES) == []
    planted, truth = rq.plant(
        rs,
        random.Random(1),
        n_date=10,
        n_cycle=5,
        n_typo=10,
        n_untag=0,
        species=SPECIES,
    )
    got = rq.qc(planted, SPECIES)
    for cls in ("date_order", "cycle", "unresolved_ref"):
        hit = {i for i, c, _ in got if c == cls}
        assert truth[cls] and truth[cls] <= hit, cls


def test_species_form_missing_from_list_is_not_a_register_error():
    rs = [rec("Cattleya Q", 2000, "Cattleya dowiana var. aurea × Cattleya mossiae")]
    assert flags_of(rs, "unresolved_ref") == set()
    assert flags_of(rs, "species_unlisted") == {0}
