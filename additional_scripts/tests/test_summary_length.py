#!usr/bin/env python3
import os
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import unittest

# test verifies that gappa_summary output:
#   [mag_id]_summary.csv
# has processed genes equal to the number of jplace
# files in the directory:
#   resources/Lenisialimosa_epa_out/{gene}/RAxML_portableTree.epa.jplace

class TestSummaryLength(unittest.TestCase):
    def test_summary_length(self):
        jplace_files = [f for f in listdir("resources/Lenisialimosa_epa_out/" + gene) if isfile(join("resources/Lenisialimosa_epa_out/" + gene, f)) and f.endswith(".jplace")]
        print(jplace_files)
        summary_file = "Lenisialimosa_summary.csv"
        gene_count = 0
        with open(summary_file, "r") as f:
            lines = []
            for line in f.readlines():
                if "Level_1" in line:
                    lines.append(line)
            for line in lines:
                count = line.split(",")[2]
                gene_count += int(count)
            print(gene_count)
            self.assertEqual(gene_count, len(jplace_files))

        
if __name__ == '__main__':
    unittest.main()