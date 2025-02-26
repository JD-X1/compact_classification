#!usr/bin/bash

echo "Working Set"
ls resources/Torulasporaglobosa_working_dataset/*fas | wc -l

echo "QFRAGS"
ls resources/Torulasporaglobosa_q_frags/*fas | wc -l

echo "MAFFT"
ls resources/Torulasporaglobosa_mafft_out/*aln | grep -v "reduced" | wc -l


echo "EPA OUT"
ls resources/Torulasporaglobosa_epa_out/*/RAxML_portableTree.epa.jplace | wc -l
