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
        # weight transitions more heavily?
        if(dict[str[i]] != dict[str[i+1]]):
            transitions[dict[str[i]]][dict[str[i + 1]]] += 9
            bases[dict[str[i]]] += 9
    init_probs[dict[str[len(str) - 1]]] += 1
    for i in range(4):
        transitions[i] /= bases[i]
    return transitions, (init_probs / len(str))


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
    num = 0
    num_5_10 = 0
    while (i < len(sequence) + 1):
        if(sequence[i - 1] == 1):
            start = i
            end = i
            while(i < len(sequence) and sequence[i] == 0):
                if(i > 4999 and i < 10001):
                    num_5_10 += 1
                else:
                    num += 1
                i += 1
                end += 1
            if(start != end):
                intervals = intervals + [(start, end)]
        i += 1
    print("conserved within nonconserved region: " +
          str(num * 100 / (len(sequence) - 5000)) + "%")
    print("conserved within conserved region: " + str(num_5_10 / 50) + "%")
    return intervals


def main():
    parser = argparse.ArgumentParser(
        description='Parse sequences into coding and noncoding regions using Phylo-HMM.')
    parser.add_argument('-seqs', action="store", dest="seqs",
                        type=str, default="data1/msa/data1-msa.fasta")
    parser.add_argument('-t', action="store", dest="transitions",
                        type=str, default="strain 17 annotations/17-bcdn.fasta")
    parser.add_argument('-m', action="store", dest="models",
                        type=str, default="hsv.nwk")
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
    #transitions = read_transitions(transition_file)
    transitions, init_probs = parse_bcdn(transition_file)
    init_models = read_models(model_file)
    init_weights = read_weights(weight_file)

    orderings = phylo.build_orderings(init_models, obs)
    phylo.reweight(orderings, init_weights)
    for _ in range(10):
        # init_probs = training.train_initial_weightings(transitions)
        emiss_probs = phylo.phylo(orderings)
        print(init_probs)
        # print(transitions)
        sequence = hmm.viterbi(transitions, emiss_probs, init_probs)
        transitions, init_probs = training.train_transitions(
            transitions, sequence)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end))
                for (start, end) in intervals]))
        f.write("\n")


if __name__ == "__main__":
    main()
