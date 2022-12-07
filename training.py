import numpy as np
import datetime


def train_transitions(transitions, sequence):
    new_transitions = np.zeros((len(transitions), len(transitions)))
    new_weightings = np.zeros(len(transitions))
    for i in range(len(sequence) - 1):
        new_transitions[sequence[i]][sequence[i+1]] += 1
        new_weightings[sequence[i]] += 1
    return new_transitions / new_weightings

# generate ancestral sequence for parent
# optimal transition matrix of child: trans_probs[x][x] = %agreements, trans_probs[x][y] = (1 - %agreements) / 3


def generate_ancestral_sequence(parent):
    # note: if both children are a -, then ancestral sequence is also a -
    
    pass


def train_model_jcm(orderings, sequence):
    # for each ordering:
    #   for each node:
    #     if it's not a leaf, set ancestral sequence
    #     set optimal branch lengths for children using ancestral sequences
    pass


def train_initial_weightings(prob_matrix):
    probs = np.zeros(len(prob_matrix))
    total = 0
    for i in range(len(prob_matrix)):
        for j in range(len(prob_matrix[0])):
            probs[j] += prob_matrix[i][j]
            total += prob_matrix[i][j]
    return (probs / total)
