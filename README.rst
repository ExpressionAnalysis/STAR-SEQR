.. image:: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR.svg?branch=master
    :target: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR

STAR-SEQR
==========
RNA Fusion Detection using the STAR-Aligner


Installation
=============
**Required**
 - biobambam2(https://github.com/gt1/biobambam2) or conda install biobambam
 - STAR(https://github.com/alexdobin/STAR) or conda install star
 - Velvet(https://github.com/dzerbino/velvet) or conda install velvet
 - samtools(https://github.com/samtools/samtools) or conda install samtools
 - UCSC utils(http://hgdownload.soe.ucsc.edu/admin/exe/) or conda install ucsc-gtftogenepred

**Install**::

    python setup.py install

Usage
======
---------------------------------
Generate a STAR index as follows:
---------------------------------
**RNA**
 - *Index*::

     STAR --runMode genomeGenerate --genomeFastaFiles hg19.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile gencodeV24lift37.gtf --runThreadN 18 --sjdbOverhang 150 --genomeSAsparseD 1

**DNA**
 - *Index*::

    STAR --runMode genomeGenerate --genomeFastaFiles hg19.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2

-------------------------
Run STAR-SEQR as follows:
-------------------------
**RNA**
 - *Align and Call*::

     starseqr.py -1 RNA_1.fastq.gz -2 RNA_2.fastq.gz -m 1 -p RNA_test -n RNA -t 12 -i path/STAR_INDEX -g gencode.gtf -r hg19.fa -vv

 - *Call Only*::

     starseqr.py -ss RNA.Chimeric.out.sam -sj RNA.Chimeric.out.junction -p RNA_test -n RNA -t 12 -i path/STAR_INDEX -g gencode.gtf -r hg19.fa -vv

**DNA**
 - *Align and Call*::

    starseqr.py -1 DNA_1.fastq.gz -2 DNA_2.fastq.gz -m 0 -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv

 - *Call Only*::

    starseqr.py DNA.Chimeric.out.sam -sj DNA.Chimeric.out.junction  -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv
