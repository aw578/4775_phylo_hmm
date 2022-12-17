import argparse
import datetime
import training
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
        strs = []
    current_str = -1
    for i in range(0, len(lines)):
        if(lines[i][0] == '>'):
            current_str += 1
        else:
            if(current_str < len(strs)):
                strs[current_str] += lines[i].strip()
            else:
                strs.append(lines[i].strip())

    nums = np.zeros((len(strs), len(strs[0])), int)
    for i in range(0, len(strs)):
        for j in range(0, len(strs[0])):
            if(strs[i][j] == 'A'):
                nums[i][j] = 0
            elif(strs[i][j] == 'C'):
                nums[i][j] = 1
            elif(strs[i][j] == 'G'):
                nums[i][j] = 2
            elif(strs[i][j] == 'T'):
                nums[i][j] = 3
            else:  # '-'
                nums[i][j] = 4

    return nums


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
        weights[i] = float(lines[i].strip())
    return weights


def parse_bcdn(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    transitions = np.zeros((4, 4))
    bases = np.zeros(4)
    init_probs = np.zeros(4)
    str = ""
    dict = {"B": 0, "C": 1, "D": 2, "N": 3}
    for i in range(0, len(lines)):
        str += lines[i].strip()
    for i in range(0, len(str) - 1):
        transitions[dict[str[i]]][dict[str[i + 1]]] += 1
        init_probs[dict[str[i]]] += 1
        bases[dict[str[i]]] += 1
    init_probs[dict[str[len(str) - 1]]] += 1
    for i in range(4):
        transitions[i] /= bases[i]

    seq = np.zeros(len(str), int)
    for i in range(0, len(str)):
        if(str[i] == 'B'):
            seq[i] = 0
        elif(str[i] == 'C'):
            seq[i] = 1
        elif(str[i] == 'D'):
            seq[i] = 2
        elif(str[i] == 'N'):
            seq[i] = 3
    return transitions, (init_probs / len(str)), seq


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


def compare_accuracy(new_seq, base_seq):
    agreements = 0
    for i in range(0, len(new_seq)):
        if(new_seq[i] == base_seq[i]):
            agreements += 1
    return (agreements / len(new_seq))


def main():
    parser = argparse.ArgumentParser(
        description='Parse sequences into coding and noncoding regions using Phylo-HMM.')
    parser.add_argument('-seqs', action="store", dest="seqs",
                        type=str, default="data2/msa/data2-msa.fasta")
    parser.add_argument('-bcdn', action="store", dest="bcdn",
                        type=str, default="strain 17 annotations/17-bcdn.fasta")
    parser.add_argument('-m', action="store", dest="models",
                        type=str, default="data2/hsv.nwk")

    args = parser.parse_args()
    fasta_file = args.seqs
    bcdn_file = args.bcdn
    model_file = args.models

    obs = read_fasta(fasta_file)
    transitions, init_probs, init_seq = parse_bcdn(bcdn_file)
    init_models = read_models(model_file)
    orderings = phylo.build_orderings(init_models, obs)
    for _ in range(2):
        emiss_probs = phylo.phylo(orderings)
        print(init_probs)
        sequence = hmm.forward_backward(transitions, emiss_probs, init_probs)
        transitions, init_probs = training.train_transitions(
            transitions, sequence)
    print(compare_accuracy(sequence, init_seq))


if __name__ == "__main__":
    main()
