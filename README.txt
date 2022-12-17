Command formatting:
python parsing.py -seqs (msa sequences) -bcdn (bcdn annotations) -m (newick trees) 

Notes on formatting:

the fasta file should be in this format:
name0
seq0
name1
seq1

the nwk file shouldbe in this format:
tree0
tree1
the names of organisms in nwk should be integers corresponding to their position in the fasta file (the organism corresponding to seq1 should be 1 in the newick tree)
use parentheses and space out commas

the transition file should be in this format:
(transitions from state 0 to states 0-n)
(transitions from state 1 to states 0-n)
...
(transitions from state n to states 0-n)

the weights file should be in this format:
weight0
weight1