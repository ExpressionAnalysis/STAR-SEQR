#!/bin/bash

STARSEQR="../starseqr.py"
cmd="python $STARSEQR -1 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_1.fastq.gz -2 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_2.fastq.gz -m 0 -p RNA_test -vv -j 1 -s 1 -d 1000000 -w 12 -i /mounts/isilon/data/indexes/DNA/h_sapiens/hg19.fa.derived/eaGene4/STAR_SEQR_hg19gencodeV24lift37_S1_RNA --ann_source gencode"
echo $cmd
eval $cmd

test_TP=$(fgrep -wf truth_rna RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$(fgrep -vwf truth_rna RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$((test_FP-1))
test_FN=$(fgrep -vwf $(cut -f2 RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt) truth | wc -l)

echo "TP = $test_TP"
echo "FP = $test_FP"
echo "FN = $test_FN"

# rm -rf RNA_test_STAR-SEQR RNA_test_STAR-SEQR.log
