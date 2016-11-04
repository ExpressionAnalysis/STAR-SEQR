.. image:: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR.svg?branch=master
    :target: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR

STAR-SEQR
==========
Description
==========
RNA Fusion Detection using the STAR-Aligner


Installation
==========
- Tools that need to be on path:
 - biobambam2
 - STAR
 - Velvet
 - Spades
 - samtools

Download a Release and install::
   
    python setup.py install

- Required python packages installed automatically:
 - intervaltree_bio
 - pandas >0.18.0
 - pysam >0.9.0
 - primer3-py


Usage
==========
Generate a STAR index as follows:
-----------
*RNA*::
     
    STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile /mounts/isilon/data/indexes/GFFs/gencodeV24lift37.gtf --runThreadN 18 --genomeSAsparseD 1

*DNA*::

    STAR --runMode genomeGenerate --genomeFastaFiles /mounts/datah/indexes/DNA/h_sapiens/hg19_scaffolds.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2

Run STAR-SEQR as follows:
-----------

*RNA-Align and Call*::

     starseqr.py -1 RNA_1.fastq.gz -2 RNA_2.fastq.gz -m 0 -p RNA_test -n RNA -t 12 -i path/STAR_INDEX --ann_source gencode -r hg19.fa -vv
 
*RNA-Call Only*::

     starseqr.py -ss RNA.Chimeric.out.sam -sj RNA.Chimeric.out.junction -p RNA_test -n RNA -t 12 -i path/STAR_INDEX --ann_source gencode -r hg19.fa -vv


*DNA-Align and Call*::

    starseqr.py -1 DNA_1.fastq.gz -2 DNA_2.fastq.gz -m 0 -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv
    
*DNA-Call Only*::

    starseqr.py DNA.Chimeric.out.sam -sj DNA.Chimeric.out.junction  -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv



