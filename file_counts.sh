#!usr/bin/bash

echo "Working Set"
ls resources/Lenisialimosa_working_dataset/*fas | wc -l

echo "QFRAGS"
ls resources/Lenisialimosa_q_frags/*fas | wc -l

echo "MAFFT"
ls resources/Lenisialimosa_mafft_out/*aln | grep -v "reduced" | wc -l


echo "EPA OUT"
ls resources/Lenisialimosa_epa_out/*/RAxML_portableTree.epa.jplace | wc -l
