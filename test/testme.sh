#!/bin/bash

STARSEQR="run_starseqr.py"
cmd="python $STARSEQR -1 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_1.fastq.gz -2 /mounts/isilon/data/eahome/q804348/opt/scripts/STAR-SEQR/test/RNA_2.fastq.gz -p RNA_test -vvv"
echo $cmd
eval $cmd

test_TP=$(fgrep -wf truth RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FP=$(fgrep -vwf truth RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt | wc -l)
test_FN=$(fgrep -vwf $(cut -f2 RNA_test_STAR-SEQR/RNA_test_STAR-SEQR_breakpoints.txt) truth | wc -l)

echo "TP = $test_TP"
echo "FP = $test_FP"
echo "FN = $test_FN"

# rm -rf RNA_test_STAR-SEQR RNA_test_STAR-SEQR.log