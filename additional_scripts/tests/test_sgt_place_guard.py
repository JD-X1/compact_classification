#!/usr/bin/env python3
"""Tests for the SGT empty/missing-jplace off-ramp (sgt_place_guard).

A per-marker placement that yields no/empty .jplace must not crash the pipeline
(`exit 2`); it is the unplaceable-marker novelty signal (dispersion_aware_placement_calling.md
sections 4-5, 7). The guard must: complete the DAG, cast zero placement votes, and
stay countable.

Run from additional_scripts/:  python -m unittest tests/test_sgt_place_guard.py
"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sgt_place_guard as g  # noqa: E402


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh)


# A minimal but schema-valid jplace (fields per Matsen et al. 2012; see
# tests/data-test/simple.jplace).
PLACED_JPLACE = {
    "tree": "((A:1{1},B:1{2})n1:1{3},C:2{4}):0{0};",
    "fields": ["edge_num", "likelihood", "like_weight_ratio",
               "distal_length", "pendant_length"],
    "placements": [{"p": [[3, -1.0, 1.0, 0.2, 0.5]], "n": ["query1"]}],
}


class TestClassifyJplace(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _path(self, name):
        return os.path.join(self.dir, name)

    def test_real_placement_is_placeable(self):
        p = self._path("ok.jplace")
        write_json(p, PLACED_JPLACE)
        status, reason = g.classify_jplace(p)
        self.assertEqual(status, "placeable")
        self.assertEqual(reason, "")

    def test_missing_file(self):
        status, reason = g.classify_jplace(self._path("does_not_exist.jplace"))
        self.assertEqual(status, "unplaceable")
        self.assertEqual(reason, "missing_file")

    def test_zero_byte_file(self):
        p = self._path("empty.jplace")
        open(p, "w").close()
        status, reason = g.classify_jplace(p)
        self.assertEqual(status, "unplaceable")
        self.assertEqual(reason, "empty_file")

    def test_malformed_json_is_unplaceable_not_crash(self):
        p = self._path("broken.jplace")
        with open(p, "w") as fh:
            fh.write("{not valid json,,,")
        status, reason = g.classify_jplace(p)
        self.assertEqual(status, "unplaceable")
        self.assertEqual(reason, "unparseable")

    def test_empty_placements_array(self):
        # The case the old `[ ! -s ]` guard let through: valid JSON, no placement.
        p = self._path("noplace.jplace")
        obj = dict(PLACED_JPLACE, placements=[])
        write_json(p, obj)
        status, reason = g.classify_jplace(p)
        self.assertEqual(status, "unplaceable")
        self.assertEqual(reason, "no_placements")

    def test_placements_key_absent(self):
        p = self._path("nokey.jplace")
        obj = {k: v for k, v in PLACED_JPLACE.items() if k != "placements"}
        write_json(p, obj)
        status, reason = g.classify_jplace(p)
        self.assertEqual(status, "unplaceable")
        self.assertEqual(reason, "no_placements")


class TestUnplaceableProfile(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.gene_dir = os.path.join(self.tmp.name, "ADK2")
        os.makedirs(self.gene_dir)
        self.profile = os.path.join(self.gene_dir, "profile.tsv")

    def tearDown(self):
        self.tmp.cleanup()

    def test_sentinel_casts_zero_votes(self):
        # "grep -v LWR" (sgt_summary) must empty the file -> no placement rows.
        g.write_unplaceable_profile(self.profile, "ADK2", "no_placements")
        self.assertTrue(os.path.exists(self.profile))
        with open(self.profile) as fh:
            data_rows = [ln for ln in fh.read().splitlines()
                         if ln.strip() and "LWR" not in ln]
        self.assertEqual(data_rows, [])

    def test_is_unplaceable_true_on_sentinel(self):
        g.write_unplaceable_profile(self.profile, "ADK2", "empty_file")
        self.assertTrue(g.is_unplaceable_profile(self.profile))

    def test_is_unplaceable_false_on_real_profile(self):
        with open(self.profile, "w") as fh:
            fh.write("LWR\tfract\taLWR\tafract\ttaxopath\n")
            fh.write("0.89\t0.90\t0.95\t0.96\tEukaryota;Amorphea;Obazoa\n")
        self.assertFalse(g.is_unplaceable_profile(self.profile))

    def test_flag_records_gene_and_reason(self):
        g.write_unplaceable_profile(self.profile, "ADK2", "no_placements")
        flag = os.path.join(self.gene_dir, "unplaceable.flag")
        self.assertTrue(os.path.exists(flag))
        with open(flag) as fh:
            content = fh.read().strip()
        self.assertEqual(content.split("\t"), ["ADK2", "no_placements"])


class TestCountUnplaceableMarkers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _make(self, gene, placed):
        d = os.path.join(self.root, gene)
        os.makedirs(d)
        p = os.path.join(d, "profile.tsv")
        if placed:
            with open(p, "w") as fh:
                fh.write("LWR\tfract\taLWR\tafract\ttaxopath\n")
                fh.write("0.9\t0.9\t0.9\t0.9\tEukaryota;X\n")
        else:
            g.write_unplaceable_profile(p, gene, "no_placements")
        return p

    def test_solarion_like_fraction(self):
        # Empirical anchor: novel Solarion ~7/232 markers unplaceable (doc section 4).
        paths = [self._make(f"g{i}", placed=(i >= 7)) for i in range(232)]
        n_unplaceable, n_total, fraction = g.count_unplaceable_markers(paths)
        self.assertEqual(n_unplaceable, 7)
        self.assertEqual(n_total, 232)
        self.assertAlmostEqual(fraction, 7 / 232, places=4)

    def test_all_placeable(self):
        paths = [self._make(f"g{i}", placed=True) for i in range(5)]
        self.assertEqual(g.count_unplaceable_markers(paths), (0, 5, 0.0))

    def test_all_unplaceable(self):
        paths = [self._make(f"g{i}", placed=False) for i in range(5)]
        self.assertEqual(g.count_unplaceable_markers(paths), (5, 5, 1.0))

    def test_empty_input_no_zero_division(self):
        self.assertEqual(g.count_unplaceable_markers([]), (0, 0, 0.0))


class TestReasonResolvedTally(unittest.TestCase):
    """The novelty-informative case (recovered + aligned + EPA ran + no placement,
    i.e. reason 'no_placements') must be separable from mere unalignable/empty
    markers (dispersion_aware_placement_calling.md sections 4-5)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _placed(self, gene):
        d = os.path.join(self.root, gene)
        os.makedirs(d)
        p = os.path.join(d, "profile.tsv")
        with open(p, "w") as fh:
            fh.write("LWR\tfract\taLWR\tafract\ttaxopath\n")
            fh.write("0.9\t0.9\t0.9\t0.9\tEukaryota;X\n")
        return p

    def _unplaceable(self, gene, reason):
        d = os.path.join(self.root, gene)
        os.makedirs(d)
        p = os.path.join(d, "profile.tsv")
        g.write_unplaceable_profile(p, gene, reason)
        return p

    def test_reason_breakdown_isolates_no_placements(self):
        paths = [self._placed(f"ok{i}") for i in range(10)]
        paths += [self._unplaceable(f"np{i}", "no_placements") for i in range(4)]
        paths += [self._unplaceable(f"ef{i}", "empty_file") for i in range(2)]
        paths += [self._unplaceable("mf0", "missing_file")]
        summary = g.summarize_unplaceable(paths)
        self.assertEqual(summary["n_markers"], 17)
        self.assertEqual(summary["n_unplaceable"], 7)
        self.assertEqual(summary["n_no_placements"], 4)
        self.assertEqual(summary["n_empty_file"], 2)
        self.assertEqual(summary["n_missing_file"], 1)
        # The novelty-informative fraction is over the no_placements flavor only.
        self.assertAlmostEqual(summary["no_placements_fraction"], 4 / 17, places=6)
        self.assertAlmostEqual(summary["unplaceable_fraction"], 7 / 17, places=6)

    def test_read_reason_from_flag(self):
        p = self._unplaceable("ADK2", "no_placements")
        self.assertEqual(g.read_unplaceable_reason(p), "no_placements")

    def test_placed_profile_has_no_reason(self):
        p = self._placed("ADK2")
        self.assertIsNone(g.read_unplaceable_reason(p))

    def test_unplaceable_without_flag_is_unknown(self):
        # Defensive: a header-only profile whose flag was lost still counts as
        # unplaceable, bucketed 'unknown', never silently dropped.
        d = os.path.join(self.root, "ORPH")
        os.makedirs(d)
        p = os.path.join(d, "profile.tsv")
        with open(p, "w") as fh:
            fh.write("LWR\tfract\taLWR\tafract\ttaxopath\n")
        summary = g.summarize_unplaceable([p])
        self.assertEqual(summary["n_unplaceable"], 1)
        self.assertEqual(summary["n_unknown"], 1)


class TestDeterminism(unittest.TestCase):
    def test_classification_and_output_reproducible(self):
        with tempfile.TemporaryDirectory() as d:
            jp = os.path.join(d, "noplace.jplace")
            write_json(jp, dict(PLACED_JPLACE, placements=[]))
            self.assertEqual(g.classify_jplace(jp), g.classify_jplace(jp))

            p1 = os.path.join(d, "a", "profile.tsv")
            p2 = os.path.join(d, "b", "profile.tsv")
            g.write_unplaceable_profile(p1, "ADK2", "no_placements")
            g.write_unplaceable_profile(p2, "ADK2", "no_placements")
            with open(p1) as fh1, open(p2) as fh2:
                self.assertEqual(fh1.read(), fh2.read())


if __name__ == "__main__":
    unittest.main()
