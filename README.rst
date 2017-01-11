|Travis| |Pypi|

=========
STAR-SEQR
=========
RNA Fusion Detection using the STAR-Aligner. Post-alignment run times are typically <10 minutes using 4 threads. DNA breakpoint detection is also supported. Development is still ongoing and several features are currently in the works.


Installation
------------

This package is tested under Linux using Python 2.7, 3.4, and 3.5.

You can install from Pypi:
::

    pip install starseqr

Or build directly from Github by cloning the project, cd into the directory and run:
::

    python setup.py install

Or from Docker:
::

    docker pull eagenomics/starseqr

Or from Bioconda:
::

    **pending approval**


**Additional Requirements**
 - biobambam2(https://github.com/gt1/biobambam2) or conda install biobambam
 - STAR(https://github.com/alexdobin/STAR) or conda install star
 - Velvet(https://github.com/dzerbino/velvet) or conda install velvet
 - samtools(https://github.com/samtools/samtools) or conda install samtools
 - UCSC utils(http://hgdownload.soe.ucsc.edu/admin/exe/) or conda install ucsc-gtftogenepred


Build a STAR Index
------------------

First make sure the dependencies are installed and generate a STAR index for your reference.

**RNA Index**
::

     STAR --runMode genomeGenerate --genomeFastaFiles hg19.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile gencodeV24lift37.gtf --runThreadN 18 --sjdbOverhang 150 --genomeSAsparseD 1

**DNA Index**
::

    STAR --runMode genomeGenerate --genomeFastaFiles hg19.fa --genomeDir ./ --runThreadN 18 --genomeSAsparseD 2


Run STAR-SEQR
---------------

STAR-SEQR can do the STAR alignments or utilize existing outputs.

RNA-Fusions
+++++++++++

*Align and Call*
::

     starseqr.py -1 RNA_1.fastq.gz -2 RNA_2.fastq.gz -m 1 -p RNA_test -n RNA -t 12 -i path/STAR_INDEX -g gencode.gtf -r hg19.fa -vv

*Call Only*::

     starseqr.py -ss RNA.Chimeric.out.sam -sj RNA.Chimeric.out.junction -p RNA_test -n RNA -t 12 -i path/STAR_INDEX -g gencode.gtf -r hg19.fa -vv

DNA-Breakpoints
+++++++++++++++

*Align and Call*
::

    starseqr.py -1 DNA_1.fastq.gz -2 DNA_2.fastq.gz -m 0 -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv

*Call Only*
::

    starseqr.py DNA.Chimeric.out.sam -sj DNA.Chimeric.out.junction  -p DNA_test -n DNA -j 2 -s 1 -t 12 -i path/STAR_INDEX_DNA --ann_source gencode -vv

Feedback
--------

Yes! Please give us your feedback, raise issues, and let us know how the tool is working for you. Pull requests are welcome.

Contributions
-------------

This project builds of the groundwork of other public contributions. Namely:

- https://github.com/pysam-developers/pysam
- https://github.com/hall-lab/svtools
- https://github.com/vishnubob/ssw
- https://github.com/libnano/primer3-py



.. |Travis| image:: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR.svg?branch=master
    :target: https://travis-ci.org/ExpressionAnalysis/STAR-SEQR

.. |Pypi| image:: https://badge.fury.io/py/starseqr.svg
    :target: https://badge.fury.io/py/starseqr

