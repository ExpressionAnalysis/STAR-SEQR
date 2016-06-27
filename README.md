# STAR-SEQR
Breakpoint and Fusion Detection

# STAR-SEQR uses STAR for alignment.Generate a genome index as follows:
# RNA
STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile /mounts/isilon/data/indexes/GFFs/gencodeV24lift37.gtf --runThreadN 18 --genomeSAsparseD 1
# DNA
STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2

# Add Paths to programs if not on path in config file.


#### Dependencies:
## Python2.7
# intervaltree_bio
# pandas # v 18.1
# pysam  # requires 0.9.0 or newer
# primer3-py

# biobambam2
# STAR
# Velvet
# Spades
# samtools



