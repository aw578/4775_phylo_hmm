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

import argparse
import numpy as np
from math import exp


def sumLogProbs(a, b):
    if a > b:
        return a + np.log1p(exp(b-a))
    else:
        return b + np.log1p(exp(a-b))


def forward_backward(trans_probs, phylo, init_probs):
    # Initialization
    # 0 = h, 1 = l
    m = len(phylo[0])
    a = len(trans_probs)
    hmm = np.zeros((a, m))

    for j in range(len(trans_probs)):
        hmm[j][0] = init_probs[j] + phylo[j][0]

    # Forwards + forwards-backwards
    for i in range(1, m):
        for k in range(0, a):
            hmm[k][i] = hmm[0][i-1] + trans_probs[0][k]
            for j in range(1, a):
                hmm[k][i] = sumLogProbs(
                    hmm[k][i], hmm[j][i-1] + trans_probs[j][k])
            hmm[k][i] += phylo[k][i]

    # Backwards + forwards-backwards
    # something about adding to hmm is fucking things up?
    back = np.zeros((a, m))
    for i in range(m - 2, -1, -1):
        for k in range(0, a):
            back[k][i] = trans_probs[k][0] + phylo[0][i + 1] + back[0][i+1]
            for j in range(1, a):
                back[k][i] = sumLogProbs(
                    back[k][i], trans_probs[k][j] + phylo[j][i + 1] + back[j][i + 1])
            hmm[k][i] += back[k][i]
        local_prob = hmm[0][i]
        for k in range(1, a):
            local_prob = sumLogProbs(local_prob, hmm[k][i])

        for k in range(0, a):
            hmm[k][i] -= local_prob
    return np.argmax(hmm, axis=0)


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


def viterbi(trans_probs, phylo, init_probs):
    # initialization
    m = len(phylo[0])
    a = len(trans_probs)
    hmm = np.zeros((a, m))
    traceback = np.zeros((a, m), int)

    # build traceback
    for j in range(0, a):
        hmm[j][0] = init_probs[j] + phylo[j][0]
    for i in range(1, m):
        for j in range(0, a):
            # find maximum:
            # TODO: np.argmax?? argmax sum
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
    sequence = np.zeros(m, int)
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
    return sequence
