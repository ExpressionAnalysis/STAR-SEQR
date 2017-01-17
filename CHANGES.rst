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
