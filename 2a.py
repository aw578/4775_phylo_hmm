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
import math
import numpy as np


'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	strs: array of stringsss
'''


def read_fasta(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        strs = ["" for _ in range(math.floor(len(lines) / 2))]
    for i in range(1, len(lines), 2):
        j = (int)((i - 1) / 2)
        strs[j] = ""
        for l in lines[i]:
            strs[j] += l.strip()
    return strs


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
    models: tree + ordering + branch length
Returns:
    likelihoods: 1 x M likelihood matrix
'''


def likelihood(model, data):
    base_conversion = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    m = len(data[0])
    likelihoods = np.zeros(m)
    for i in range(0, m):
        # temp
        likelihoods[i] = model[base_conversion[data[0][i]]]
    return likelihoods


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    models: array of A models (trees + orderings + branch lengths)
    data: array of N M-char sequences
Returns:
    likelihoods: A x M log likelihood matrix
'''


def phylo(models, data):
    a = len(models)
    m = len(data[0])
    likelihoods = np.zeros((a, m))
    for i in range(0, a):
        likelihoods[i] = likelihood(models[i], data)
    return likelihoods


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
    n = len(obs)
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


''' Returns a list of non-overlapping intervals describing the GC rich regions.
Arguments:
	sequence: list of hidden states
Returns:
	intervals: list of tuples (i, j), 1 <= i <= j <= len(sequence), that
                describe GC rich regions in the input list of hidden states.
'''


def find_intervals(sequence):
    intervals = []
    i = 1
    while (i < len(sequence) + 1):
        if(sequence[i - 1] == 0):
            start = i
            end = i
            while(i < len(sequence) and sequence[i] == 0):
                i += 1
                end += 1
            intervals = intervals + [(start, end)]
        i += 1
    return intervals


def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Viterbi.')
    parser.add_argument('-f', action="store", dest="f",
                        type=str, default="hmm-sequence.fasta")
    parser.add_argument('-mu', action="store", dest="mu",
                        type=float, default=0.01)
    parser.add_argument('-out', action="store", dest="out",
                        type=str, default="viterbi-intervals.txt")

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu
    intervals_file = args.out

    obs_sequence = read_fasta(fasta_file)
    transition_probabilities = np.array([
        [np.log(1 - mu), np.log(mu)],
        [np.log(mu), np.log(1 - mu)]])
    emission_probabilities = np.array([
        [np.log(0.13), np.log(0.37), np.log(0.37), np.log(0.13)],
        [np.log(0.32), np.log(0.18), np.log(0.18), np.log(0.32)]
    ])
    emission_probabilities = phylo(emission_probabilities, obs_sequence)
    initial_probabilities = np.array([np.log(0.5), np.log(0.5)])
    sequence, p = viterbi(obs_sequence, transition_probabilities,
                          emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")
    print("Viterbi probability: {:.2f}".format(p))


if __name__ == "__main__":
    main()
