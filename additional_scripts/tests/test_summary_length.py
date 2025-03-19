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
        jplace_files = listdir("resources/Lenisialimosa_epa_out/")
        print("Jplace files: ")
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
    def test_pipeline_geneCount_check(self):
        working_dataset_list = listdir("resources/Lenisialimosa_working_dataset/")
        spitter_dataset_list = listdir("resources/Lenisialimosa_q_frags/")
        mafft_dataset_list = listdir("resources/Lenisialimosa_mafft_out/")
        self.assertEqual(len(working_dataset_list), len(spitter_data_list))
        self.assertEqual(len(splitter_dataset_list), len(mafft_data_list))

        
if __name__ == '__main__':
    unittest.main()