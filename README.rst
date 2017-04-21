|Travis| |Pypi| |Conda|

=========
STAR-SEQR
=========
RNA Fusion Detection and Quantification using STAR.

Post-alignment run times are typically <20 minutes using 4 threads.  Development is still ongoing and several features are currently in the works. DNA breakpoint detection is still experimental.


Installation
------------

This package is tested under Linux using Python 2.7, 3.4, and 3.5.

You can install from Pypi. Please use a recent version of pip and cython:
::

    pip install -U pip
    pip install -U cython
    pip install starseqr

Or build directly from Github by cloning the project, cd into the directory and run:
::

    python setup.py install

Or from Docker:
::

    docker pull eagenomics/starseqr

Or from Bioconda:
::

    conda install starseqr


**Additional Requirements**
 - biobambam2(https://github.com/gt1/biobambam2) or conda install biobambam
 - STAR(https://github.com/alexdobin/STAR). Must use >2.5.3a. conda install star
 - Velvet(https://github.com/dzerbino/velvet) or conda install velvet
 - samtools(https://github.com/samtools/samtools) or conda install samtools
 - Salmon(https://combine-lab.github.io/salmon/) or conda install salmon
 - UCSC utils(http://hgdownload.soe.ucsc.edu/admin/exe/) or conda install ucsc-gtftogenepred
 - gffread(http://ccb.jhu.edu/software/stringtie/dl/gffread-0.9.8c.tar.gz) or conda install gffread


Build a STAR Index
------------------

First make sure the dependencies are installed and generate a STAR index for your reference.

**RNA Index**
::

     STAR --runMode genomeGenerate --genomeFastaFiles hg19.fa --genomeDir STAR_SEQR_hg19gencodeV24lift37_S1_RNA --sjdbGTFfile gencodeV24lift37.gtf --runThreadN 18 --sjdbOverhang 150 --genomeSAsparseD 1


Run STAR-SEQR
--------------

STAR-SEQR can perform alignment or utilize existing outputs from STAR. Note- STAR-SEQR alignment parameters have been tuned for fusion calling.


**Align and Call**
::

     starseqr.py -1 RNA_1.fastq.gz -2 RNA_2.fastq.gz -m 1 -p RNA_test -n RNA -t 12 -i path/STAR_INDEX -g gencode.gtf -r hg19.fa -vv


Outputs
-------

Breakpoints.txt and Candidates.txt have the following columns:

+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Values**          | **Description**                                                                                                                                                        |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NAME                | Gene Symbols for left and right fusion partners                                                                                                                        |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NREAD_SPANS         | The number of paired reads that are discordant spanning and suppor the fusion                                                                                          |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NREAD_JXNLEFT       | The number of paired reads that are anchored on the left side of the gene fusion                                                                                       |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NREAD_JXNRIGHT      | The number of paired reads that are anchored on the right side of the gene fusion                                                                                      |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| FUSION_CLASS        | Classification of fusion based on chromosomal location, distance and strand. [GENE_INTERNAL, TRANSLOCATION, READ_THROUGH, INTERCHROM_INVERTED, INTERCHROM_INTERSTRAND] |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| SPLICE_TYPE         | Classification of the fusion breakpoint. If on the exon boundary is CANONICAL, else NON-CANONICAL                                                                      |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| BRKPT_LEFT          | The 0-based genomic position of the fusion breakpoint for the left gene partner                                                                                        |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| BRKPT_RIGHT         | The 0-based genomic position of the fusion breakpoint for the right gene partner                                                                                       |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| LEFT_SYMBOL         | The left gene symbol                                                                                                                                                   |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| RIGHT_SYMBOL        | The right gene symbol                                                                                                                                                  |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ANNOT_FORMAT        | The description of keys that are used in the ANNOT column. Similar to VCF FORMAT notation.                                                                             |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| LEFT_ANNOT          | The values described in the ANNOT_FORMAT column for the left gene breakpoint                                                                                           |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| RIGHT_ANNOT         | The values described in the ANNOT_FORMAT column for the right gene breakpoint                                                                                          |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| DISTANCE            | The genomic distance between breakpoints. Empty if a translocation.                                                                                                    |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ASSEMBLED_CONTIGS   | The velvet assembly of the supporting chimeric reads                                                                                                                   |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ASSEMBLY_CROSS_JXN  | A boolean value indicating if the assembly crosses the putative breakpoint                                                                                             |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| PRIMERS             | Primers left, right designed against the highest expressing predicted fusion transcript                                                                                |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ID                  | Internal notation of STAR-SEQR breakpoints.                                                                                                                            |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| SPAN_CROSSHOM_SCORE | Homology score with range of [0-1] to indicate the probability of spanning chimeric reads mapping to both gene partners                                                |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| JXN_CROSSHOM_SCORE  | Homology score with range of [0-1] to indicate the probability of junction chimeric reads mapping to both gene partners                                                |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| OVERHANG_DIVERSITY  | The number of unique fragments that fall from left anchored split-reads onto the right gene and vice-versa.                                                            |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MINFRAG20           | The number of overhang fragments that have at least 20 bases                                                                                                           |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MINFRAG35           | The number of overhang fragments that have at least 35 bases                                                                                                           |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TPM_FUSION          | Expression of the most abundant fusion transcript expressed in transcripts per million                                                                                 |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TPM_LEFT            | Expression of the most abundant left transcript expressed in transcripts per million                                                                                   |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TPM_RIGHT           | Expression of the most abundant right transcript expressed in transcripts per million                                                                                  |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MAX_TRX_FUSION      | Highest expressing fusion transcript. Expression corresponds to TPM_FUSION                                                                                             |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| DISPOSITION         | Values to indicate PASS or other specific reasons for failure                                                                                                          |
+---------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

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

.. |Conda| image:: https://anaconda.org/bioconda/starseqr/badges/installer/conda.svg
    :target: https://bioconda.github.io/recipes/starseqr/README.html
