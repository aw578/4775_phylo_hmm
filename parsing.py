import argparse
import math
import numpy as np
import hmm
import phylo

'''Reads the fasta file and outputs the sequences to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	strs: array of sequences
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


'''
Arguments:
  filename: the file to read
Return:
  transitions: A x A matrix of transition probabilities
'''


def read_transitions(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        transitions = np.zeros((len(lines), len(lines)))
    for i in range(0, len(lines)):
        transitions[i] = np.array(lines[i].split())
    return transitions


'''
Arguments:
  filename: the file to read
Return:
  weights: 1 x A matrix of weightings
'''


def read_weights(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        weights = np.zeros(len(lines))
    for i in range(0, len(lines)):
        weights[i] = int(lines[i].strip())
    return weights


'''
Arguments:
  filename: the file to read
Return:
  trees: 1 x A list of newick trees
'''


def read_models(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        trees = []
    for i in range(0, len(lines)):
        trees.append(lines[i].strip())
    return trees


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


def find_init_probs(prob_matrix):
    probs = np.zeros(len(prob_matrix))
    total = 0
    for i in range(len(prob_matrix)):
        for j in range(len(prob_matrix[0])):
            probs[j] += prob_matrix[i][j]
            total += prob_matrix[i][j]
    return (probs / total)


def main():
    parser = argparse.ArgumentParser(
        description='Parse sequences into coding and noncoding regions using Phylo-HMM.')
    parser.add_argument('-seqs', action="store", dest="seqs",
                        type=str, default="hmm-sequence.fasta")
    parser.add_argument('-t', action="store", dest="transitions",
                        type=str, default="transitions.txt")
    parser.add_argument('-m', action="store", dest="models",
                        type=str, default="models.nwk")
    parser.add_argument('-w', action="store", dest="weights",
                        type=str, default="weights.txt")
    parser.add_argument('-out', action="store", dest="out",
                        type=str, default="viterbi-intervals.txt")

    args = parser.parse_args()
    fasta_file = args.seqs
    transition_file = args.transitions
    model_file = args.models
    weight_file = args.weights
    intervals_file = args.out

    obs = read_fasta(fasta_file)
    transitions = read_transitions(transition_file)
    init_models = read_models(model_file)
    init_weights = read_weights(weight_file)

    init_probs = find_init_probs(transitions)
    orderings = phylo.build_orderings(init_models)
    phylo.reweight(orderings, init_weights)

    emiss_probs = phylo.phylo(orderings, obs)

    sequence, p = hmm.viterbi(obs, transitions, emiss_probs, init_probs)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")
    print("Viterbi probability: {:.2f}".format(p))


if __name__ == "__main__":
    main()
