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
