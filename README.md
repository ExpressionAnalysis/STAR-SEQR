# STAR-SEQR
Breakpoint and Fusion Detection. More description to follow.

### Generate STAR index as follows:
* RNA
  * STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile /mounts/isilon/data/indexes/GFFs/gencodeV24lift37.gtf --runThreadN 18 --genomeSAsparseD 1

* DNA
  * STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2


### Dependencies:
* Python2.7
  * intervaltree_bio
  * pandas # v 18.1
  * pysam  # requires 0.9.0 or newer

* Other tools
  * primer3-py
  * biobambam2
  * STAR
  * Velvet
  * Spades
  * samtools

### Config file must have paths to programs


