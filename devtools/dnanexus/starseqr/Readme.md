<!-- dx-header -->
# STAR-SEQR (DNAnexus Platform App)

RNA Fusion Detection and Quantification

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

<!-- Insert a description of your app here -->
https://github.com/ExpressionAnalysis/STAR-SEQR

Two run modes:
1) Do alignment and fusion calling. Specify inputs:
   * fastq1
   * fastq2
   * star_index
   * genome_fasta
   * transcript_gtf
   * prefix

2) Use existing alignment to do fusion calling. Specify inputs:
   * fastq1
   * fastq2
   * genome_fasta
   * transcript_gtf
   * chimeric_juncs
   * chimeric_sam | chimeric_bam
   * prefix
