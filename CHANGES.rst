STAR-SEQR 0.1.1:
 * STAR-SEQR is now compatible with Py3
 * Travis CI for python versions 2.7, 3.4, 3.5
 * Checks to ensure primer design targets are compatible w/seq length
 * Annotation field is now ordered by junction dist, cdslen and transcript name for consistency

STAR-SEQR 0.2.0:
 * The smallest portion of a split read is now used to calculate number of frags >20 as MINFRAG20 or MINFRAG35
 * More filters are applied to remove FPs with certain read characteristics such as diversity, MINFRAG20, and Cross-Homology
 * Candidate info has been converted to not contain nested lists
 * New output format with more info

STAR-SEQR 0.2.1:
 * Added a script for making a DREAM challenge compliant bedpe

STAR-SEQRv0.2.2:
 * Change absolute path requirements for gtf and fasta to accomodate docker

STAR-SEQRv0.2.3:
 * Changed all realpath statements to abspath.

STAR-SEQRv0.2.4:
 * Updated to work with pandas 19.2 and numpy 11.3
 * Fasta index is created if non-existent
 * Candidates.txt shows PASS or fail reasons

STAR-SEQRv0.3.0:
 * Big release!
 * Better modularity in code
 * Speed improvements for all annotation functions
 * Improved logging
 * Fixed a bug in bed subsetting when using the "both" style
 * No longer produces the same style of bedpe/VCF for Fusions.
 * New Fusion bedpe using proper coordinates
 * Support fastqs are broken into more coherent groupings: span, split, overhang
 * Cross homology checks now use span and overhang reads separately
 * Fixed a primer bug where breakpoints were not being used correctly in some cases
 * Fixed several bugs that were introduced that prevented DNA mode from working.

STAR-SEQRv0.3.1:
 * Fixed install script and now use sdist and twine for Pypi uploads.

STAR-SEQRv0.4.0:
 * Updates to STAR functions to use the new chimMainSegmentMultNmax parameter.
 * Allow more ways to run STAR-SEQR from different STAR outputs.
 * Salmon is now used to quantify candidate fusions and partner transcripts.
 * Additional parameter (-x) to pass in a transcript reference for more stable counts.
 * Greater sensitivity. Now rescues reads with 1 read support as long as all other filters and expression metrics are met.
 * Uses more filtering approaches including basequalities, expression of transcripts relative to fusion,
 * Fixed a bug so fusion junctions are now normalized to 0-base genomic coordinates.

STAR-SEQRv0.5.0:
 * Creates a central folder "chimeric_transcripts" for all fasta files.
 * Functions now pull fastas from central chimeric transcript folder
 * Reference transcripts are now extracted directly from the GTF. The (-x) parameter is no longer in use.
 * Multimapping homologous fusions are considered through a network graph approach where homology and expression are used to prune clusters. If homologous the highest expressing homolog passes the filter. If they share a breakpoint but are not considered homologoues, then all partners pass the filter.
 * The highest expressing fusion transcript is now part of the output
 * Primers are now generated from the highest expressing transcript.

STAR-SEQRv0.5.1:
 * Added bioconda recipe
 * Updated basespace script
 * Modified threshold for num of jxns to minfrag20 ratio from 10% to 1%. This filter was removing TPs in some datasets and other filters are catching FPs in our training datasets.

 STAR-SEQRv0.6.0:
 * With a compendia of training data tweaked thresholds to keep FP from exploding in especially noisy samples such as FFPE.
 * Added more info columns to the output

STAR-SEQRv0.6.1:
 * Fixed a bug where fastq input was marked required but would break when run using existing alignments
 * Added cwl for starseqr

STAR-SEQRv0.6.2:
 * Changed map multiprocessing to map_async which leads to overall better performance
 * Filters mitochondria junctions earlier in the process which also imporoves performance
 * symlink gtf into folder and make table locally rather than externally. Helpful for CWL
 * changed memory to low requirement for pandas read_csv to avoid segfaults on large inputs
 * use salmon "A" library type

STAR-SEQRv0.6.3:
 * fixed a filtering step that would die if all junctions were based on previous filtering

STAR-SEQRv0.6.4:
 * Allow prefix to be a string or a full path. Helpful for writing outside of current working directory.
 * Removed DNA functions as this has been deprecated for some time.

STAR-SEQRv0.6.5:
 * Primer generation step was taking too long. Modified to search only 200bp on each side of breakpoint.

STAR-SEQRv0.6.6:
 * If no fusion transcript is identified skip primer generation. This was causing errors in some systems.

