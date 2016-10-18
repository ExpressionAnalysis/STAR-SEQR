# STAR-SEQR
Breakpoint and Fusion Detection.


#Installation
* Tools that need to be on path
  * biobambam2
  * STAR
  * Velvet
  * Spades
  * samtools

* STAR-SEQR
  * Download a release
  * cd into the directory
  * python setup.py build
  * python setup.py install
  * python setup.py clean

### Dependencies:
* Python2.7
  * intervaltree_bio
  * pandas >0.18.0
  * pysam >0.9.0
  * primer3-py

#Usage
### Generate STAR index as follows:
* RNA
  * STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile /mounts/isilon/data/indexes/GFFs/gencodeV24lift37.gtf --runThreadN 18 --genomeSAsparseD 1

* DNA
  * STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2

### Run STAR-SEQR as follows:
* RNA
  * Align and Call:
    * starseqr.py -1 RNA_1.fastq.gz -2 RNA_2.fastq.gz -m 0 -p RNA_test -n RNA -t 12 -i path/STAR_INDEX --ann_source gencode -r hg19.fa -vv
  * Call Only:
    * starseqr.py -ss RNA.Chimeric.out.sam -sj RNA.Chimeric.out.junction -p RNA_test -n RNA -t 12 -i path/STAR_INDEX --ann_source gencode -r hg19.fa -vv

* DNA
  * Align and Call:
    * starseqr.py -1 DNA_1.fastq.gz -2 DNA_2.fastq.gz -m 0 -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv
  * Call Only:
    * starseqr.py DNA.Chimeric.out.sam -sj DNA.Chimeric.out.junction  -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv



