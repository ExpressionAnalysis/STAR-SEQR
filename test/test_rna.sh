#!/bin/bash

STARSEQR="../starseqr.py"
cmd="python $STARSEQR -1 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_1.fastq.gz -2 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_2.fastq.gz -m 1 -p RNA_test -t 12 -i /mounts/isilon/data/indexes/DNA/h_sapiens/hg19.fa.derived/eaGene4/STAR_SEQR_hg19gencodeV24lift37_S1_RNA --ann_source gencode -r /ea/indexes/DNA/h_sapiens/hg19.fa -vv"
echo $cmd
eval $cmd

test_TP=$(fgrep -wf truth_rna RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$(fgrep -vwf truth_rna RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$((test_FP-1))
cmd2="cut -f15 RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt > tmp_rna"
eval $cmd2
test_FN=$(fgrep -vwf tmp_rna truth_rna | wc -l)

echo "TP = $test_TP"
echo "FP = $test_FP"
echo "FN = $test_FN"

# rm -rf RNA_test_STAR-SEQR RNA_test_STAR-SEQR.log
