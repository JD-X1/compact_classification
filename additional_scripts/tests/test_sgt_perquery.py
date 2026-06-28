#!/usr/bin/env python3
"""TDD suite for the SGT->taxonomy_report adapter (sgt_perquery).

In SGT mode each marker gene is one query; its best-hit taxopath is that marker's
vote. taxonomy_report.py computes the between-marker dispersion D1 (effective_num_taxa;
Hill q=1 number, Hill 1973 / Jost 2006) but needs a per-query table (name+LWR+taxopath)
and errors on a bare profile.tsv. This adapter reshapes the per-gene best-hit profiles
into that table so D1 is produced for SGT, WITHOUT changing gappa flags or the concat path.

Run from additional_scripts/:  python -m pytest tests/test_sgt_perquery.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sgt_perquery as A          # noqa: E402
import sgt_place_guard as G       # noqa: E402
import taxonomy_report as TR      # noqa: E402

HEADER = "LWR\tfract\taLWR\tafract\ttaxopath"


def _profile(gene_dir, lwr, taxopath, header=HEADER):
    os.makedirs(gene_dir, exist_ok=True)
    p = os.path.join(gene_dir, "profile.tsv")
    cols = header.split("\t")
    row = {"LWR": str(lwr), "fract": "0", "aLWR": "0", "afract": "0", "taxopath": taxopath}
    with open(p, "w") as fh:
        fh.write(header + "\n")
        fh.write("\t".join(row[c] for c in cols) + "\n")
    return p


# ---- A: vote assembly -----------------------------------------------------
def test_A1_placeable_gene_yields_one_record(tmp_path):
    p = _profile(str(tmp_path / "ADK2"), 0.9, "Eukaryota;Chlorophyta")
    rows = A.assemble_per_query([p])
    assert rows == [("ADK2", "0.9", "Eukaryota;Chlorophyta")]

def test_A2_unplaceable_sentinel_contributes_nothing(tmp_path):
    d = str(tmp_path / "AGB1")
    G.write_unplaceable_profile(os.path.join(d, "profile.tsv"), "AGB1", "no_placements")
    assert A.assemble_per_query([os.path.join(d, "profile.tsv")]) == []

def test_A3_columns_resolved_by_header_not_position(tmp_path):
    # taxopath before LWR: adapter must key off header names.
    hdr = "taxopath\tfract\taLWR\tafract\tLWR"
    p = _profile(str(tmp_path / "MCM-B"), 0.77, "Eukaryota;Ciliophora", header=hdr)
    rows = A.assemble_per_query([p])
    assert rows == [("MCM-B", "0.77", "Eukaryota;Ciliophora")]

def test_A4_gene_taken_from_parent_dir(tmp_path):
    p1 = _profile(str(tmp_path / "g1"), 0.5, "Eukaryota;A")
    p2 = _profile(str(tmp_path / "g2"), 0.6, "Eukaryota;B")
    genes = {r[0] for r in A.assemble_per_query([p1, p2])}
    assert genes == {"g1", "g2"}


# ---- B: output is taxonomy_report-consumable ------------------------------
def test_B5_written_table_feeds_read_profile(tmp_path):
    p1 = _profile(str(tmp_path / "g1"), 0.9, "Eukaryota;Chlorophyta")
    p2 = _profile(str(tmp_path / "g2"), 0.8, "Eukaryota;Ciliophora")
    out = str(tmp_path / "sgt_per_query.tsv")
    A.write_per_query_tsv([p1, p2], out)
    gappa = TR.read_profile(out)              # must NOT raise pick_query_col error
    assert set(gappa) == {"g1", "g2"}
    assert gappa["g1"]["lineage"] == "Eukaryota;Chlorophyta"
    assert abs(gappa["g1"]["lwr"] - 0.9) < 1e-9


# ---- C: D1 correctness (statistical validation) ---------------------------
def _heterorows(tmp_path, lineages, min_markers=10):
    profiles = [_profile(str(tmp_path / f"g{i}"), 1.0, lin) for i, lin in enumerate(lineages)]
    out = str(tmp_path / "pq.tsv")
    A.write_per_query_tsv(profiles, out)
    gappa = TR.read_profile(out)
    return TR.heterogeneity_rows(gappa, "mag", min_markers=min_markers)

def test_C6_consensus_effective_near_one(tmp_path):
    rows = _heterorows(tmp_path, ["Eukaryota;Chlorophyta"] * 12)
    deepest = [r for r in rows if r["level"] == 1][0]
    assert deepest["richness"] == 1
    assert abs(deepest["effective_num_taxa"] - 1.0) < 1e-3

def test_C7_even_split_effective_near_two(tmp_path):
    rows = _heterorows(tmp_path,
                       ["Eukaryota;Chlorophyta"] * 10 + ["Eukaryota;Ciliophora"] * 10)
    top = [r for r in rows if r["level"] == 1][0]
    assert top["richness"] == 2
    assert 1.9 <= top["effective_num_taxa"] <= 2.0


# ---- D: off-ramp interplay -------------------------------------------------
def test_D8_unplaceable_excluded_from_placed_denominator(tmp_path):
    placeable = [_profile(str(tmp_path / f"g{i}"), 1.0, "Eukaryota;Chlorophyta")
                 for i in range(10)]
    unplaceable = []
    for i in range(5):
        d = str(tmp_path / f"u{i}")
        G.write_unplaceable_profile(os.path.join(d, "profile.tsv"), f"u{i}", "no_placements")
        unplaceable.append(os.path.join(d, "profile.tsv"))
    out = str(tmp_path / "pq.tsv")
    A.write_per_query_tsv(placeable + unplaceable, out)
    gappa = TR.read_profile(out)
    rows = TR.heterogeneity_rows(gappa, "mag", min_markers=10)
    assert rows[0]["n_markers_placed"] == 10        # 5 unplaceable did not vote
