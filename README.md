This script is the first in a set of custom scripts to facilitate in-silico screening of 
sequences for exact PCR primers and probe hybridization sites to screen for strains of target
organism that may not be detected by lab-developed PCR assays. 
This script searches a nucleotide fasta file for an exact match to the PCR amplicon, using information
contained in an assay-specific oligonucleotide (csv) definition file. The amplicon is created
from primers and probe in the correct order, with intervening regions of 0-100 N's.
Each query fasta sequence is read as a Biopython SeqRecord and searched for an exact match to the amplicon
(or its reverse complement) - allowing degenerate bases - and sequences lacking an exact match are output
to fasta for downstream analysis. A summary containing record IDs of non-matching sequences is printed to a text file.
