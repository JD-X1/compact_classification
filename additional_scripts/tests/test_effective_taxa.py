#!/usr/bin/env python3
"""Tests for the effective-number-of-taxa metric (Upgrade 1) in taxonomy_report.

Run from additional_scripts/:  python -m unittest tests/test_effective_taxa.py
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import taxonomy_report as tr  # noqa: E402

BIG = 1e9  # n_eff large enough that the Miller-Madow term vanishes


class TestShannonEffective(unittest.TestCase):
    def test_single_lineage_is_one(self):
        H, H_mm, richness, eff = tr.shannon_effective({"Chlorophyta": 1.0}, BIG)
        self.assertEqual(richness, 1)
        self.assertAlmostEqual(H, 0.0, places=9)
        self.assertAlmostEqual(eff, 1.0, places=6)

    def test_even_binary_split_is_two(self):
        H, H_mm, richness, eff = tr.shannon_effective({"A": 0.5, "B": 0.5}, BIG)
        self.assertEqual(richness, 2)
        self.assertAlmostEqual(H, math.log(2), places=9)
        self.assertAlmostEqual(eff, 2.0, places=6)

    def test_lopsided_split_collapses_toward_one(self):
        # 28 vs 2 of equal weight -> p = {0.9333, 0.0667}
        support = {"major": 28 / 30, "minor": 2 / 30}
        H, H_mm, richness, eff = tr.shannon_effective(support, BIG)
        self.assertEqual(richness, 2)
        self.assertAlmostEqual(eff, math.exp(H), places=6)
        self.assertTrue(1.0 < eff < 1.4, f"expected ~1.28, got {eff}")

    def test_effective_capped_at_richness(self):
        # Tiny n_eff inflates H_mm; effective number must not exceed richness.
        H, H_mm, richness, eff = tr.shannon_effective({"A": 0.5, "B": 0.5}, n_eff=2.0)
        self.assertGreater(H_mm, H)
        self.assertEqual(eff, 2.0)


class TestWeightedSupport(unittest.TestCase):
    def test_taxon_at_depth(self):
        self.assertEqual(tr.taxon_at_depth("Euk;Chlorophyta;Chlorophyceae", 1), "Chlorophyta")
        self.assertEqual(tr.taxon_at_depth("Euk;Chlorophyta", 5), "")

    def test_weighting_by_confidence(self):
        # Two markers on taxon A (weights 0.9, 0.9), one on B (weight 0.1).
        obs = [("A", 0.9), ("A", 0.9), ("B", 0.1)]
        support, n, n_eff = tr.weighted_support_at_depth(obs, 0)
        self.assertEqual(n, 3)
        self.assertAlmostEqual(support["A"], 1.8 / 1.9, places=6)
        self.assertAlmostEqual(support["B"], 0.1 / 1.9, places=6)
        # Kish n_eff = (1.9)^2 / (0.81+0.81+0.01) = 3.61 / 1.63
        self.assertAlmostEqual(n_eff, (1.9 ** 2) / 1.63, places=6)

    def test_zero_weight_markers_ignored(self):
        obs = [("A", 0.0), ("B", 0.5)]
        support, n, n_eff = tr.weighted_support_at_depth(obs, 0)
        self.assertEqual(n, 1)
        self.assertEqual(set(support), {"B"})


class TestHeterogeneityRows(unittest.TestCase):
    def _gappa(self, items):
        return {f"MAG_{i}": {"lineage": lin, "lwr": w}
                for i, (lin, w) in enumerate(items)}

    def test_clean_mag_effective_near_one(self):
        gappa = self._gappa([("Euk;Chlorophyta", 1.0)] * 20)
        rows = tr.heterogeneity_rows(gappa, "cleanMAG", min_markers=10)
        deepest = [r for r in rows if r["level"] == 1][0]
        self.assertEqual(deepest["status"], "ok")
        self.assertEqual(deepest["richness"], 1)
        self.assertAlmostEqual(deepest["effective_num_taxa"], 1.0, places=4)

    def test_binary_chimera_effective_near_two(self):
        gappa = self._gappa(
            [("Euk;Chlorophyta", 1.0)] * 10 + [("Bac;Proteobacteria", 1.0)] * 10
        )
        rows = tr.heterogeneity_rows(gappa, "chimera", min_markers=10)
        top = [r for r in rows if r["level"] == 0][0]
        self.assertEqual(top["richness"], 2)
        self.assertTrue(1.9 <= top["effective_num_taxa"] <= 2.0)

    def test_low_marker_count_flagged(self):
        gappa = self._gappa([("Euk;Chlorophyta", 1.0)] * 3)
        rows = tr.heterogeneity_rows(gappa, "tiny", min_markers=10)
        self.assertTrue(all(r["status"] == "insufficient_markers" for r in rows))

    def test_unplaced_markers_excluded(self):
        gappa = self._gappa([("Euk;Chlorophyta", 1.0), ("", 0.0), ("  ", 0.5)])
        rows = tr.heterogeneity_rows(gappa, "m", min_markers=1)
        self.assertEqual(rows[0]["n_markers_placed"], 1)


if __name__ == "__main__":
    unittest.main()
