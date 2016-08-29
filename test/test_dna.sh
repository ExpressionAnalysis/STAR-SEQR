#!/bin/bash

STARSEQR="../starseqr.py"
cmd="python ../starseqr.py -1 DNA_1.fastq.gz -2 DNA_2.fastq.gz -m 0 -p DNA_test -n DNA -vv -j 2 -s 2 -d 1000000 -w 12 -i /mounts/isilon/data/indexes/DNA/h_sapiens/hg19.fa.derived/eaGene4/STAR_SEQR_HUMAN_S2 --ann_source gencode"
echo $cmd
eval $cmd


test_TP=$(fgrep -wf truth_dna DNA_test_STAR-SEQR/DNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$(fgrep -vwf truth_dna DNA_test_STAR-SEQR/DNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$((test_FP-1))


echo "TP = $test_TP"
echo "FP = $test_FP"

