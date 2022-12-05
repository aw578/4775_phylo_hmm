#!/usr/bin/env python3

'''Script for computing GC-rich and GC-poor intervals in a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states
    -out: file to output intervals to (1 interval per line)

Outputs:
    File with list of intervals (a_i, b_i) such that bases a_i to b_i are
    classified as GC-rich.
    
Example Usage:
    python 2a.py -f hmm-sequence.fasta -mu 0.01 -out viterbi-intervals.txt
'''

import datetime
import numpy as np


''' Outputs the Viterbi decoding of a given observation.
Arguments:
	obs: list of N M-char sequences
	trans_probs: A x A matrix of transition log probs
	phylo: A x M log likelihood matrix
	init_probs: A x 1 matrix of initial log probs
Returns:
	l: list of most likely hidden states at each position
        (list of hidden states)
	p: log-probability of the returned hidden state sequence
'''


def viterbi(obs, trans_probs, phylo, init_probs):
    # initialization
    m = len(obs[0])
    a = len(trans_probs)
    hmm = np.zeros((a, m))
    traceback = np.zeros((a, m), int)

    # build traceback
    for j in range(0, a):
        hmm[j][0] = init_probs[j] + phylo[j][0]
    for i in range(1, m):
        for j in range(0, a):
            # find maximum:
            max_state = 0
            max_state_prob = hmm[0][i - 1] + trans_probs[0][j]
            for k in range(1, a):
                new_state_prob = hmm[k][i - 1] + trans_probs[k][j]
                if(new_state_prob > max_state_prob):
                    max_state = k
                    max_state_prob = new_state_prob
            hmm[j][i] = max_state_prob + phylo[j][i]
            traceback[j][i] = max_state
    # find starting point
    sequence = np.zeros(m)
    curr_state = 0
    curr_state_prob = hmm[0][m - 1]
    for k in range(1, a):
        new_state_prob = hmm[k][m - 1]
        if(new_state_prob > curr_state_prob):
            curr_state = k
            curr_state_prob = new_state_prob

    sequence[m - 1] = curr_state

    # reconstruct maximum from traceback
    for i in range(m - 1, 0, -1):
        curr_state = traceback[curr_state][i]
        sequence[i - 1] = curr_state
    return sequence, curr_state_prob
