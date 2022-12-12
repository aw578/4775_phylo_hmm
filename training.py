import numpy as np
import datetime


def train_transitions(transitions, sequence):
    new_transitions = np.zeros((len(transitions), len(transitions)))
    new_weightings = np.zeros(len(transitions))
    new_probs = np.zeros(len(transitions))
    for i in range(len(sequence) - 1):
        new_transitions[sequence[i]][sequence[i+1]] += 1
        new_weightings[sequence[i]] += 1
        new_probs[sequence[i]] += 1
    new_probs[sequence[len(sequence) - 1]] += 1
    for i in range(len(new_transitions)):
        new_transitions[i] = new_transitions[i] / new_weightings[i]
    return new_transitions, (new_probs / len(sequence))

# generate ancestral sequence for parent
# optimal transition matrix of child: trans_probs[x][x] = %agreements, trans_probs[x][y] = (1 - %agreements) / 3


def set_ancestral_sequence(parent):
    ancestral_sequence = np.argmax(parent.probs, axis=1)
    # if both children are -, then ancestral sequence is also a -
    for i in range(0):
        if(parent.left.ancestral_sequence[i] == '-' and parent.right.ancestral_sequence[i] == '-'):
            ancestral_sequence[i] = 4

    parent.ancestral_sequence = ancestral_sequence


def set_jcm_branch_lengths(parent):
    # percentage of agreements between ancestral sequence and child sequence on diagonals
    agreements = 0
    p_seq = parent.ancestral_sequence
    l_seq = parent.left.ancestral_sequence
    r_seq = parent.right.ancestral_sequence
    for i in range(len(parent.ancestral_sequence)):
        agreements += (p_seq[i] == '-' or l_seq[i] ==
                       '-' or p_seq[i] == l_seq[i])
    agreements /= len(parent.ancestral_sequence)
    new_jcm = np.full((4, 4), (1 - agreements) / 3)
    np.fill_diagonal(new_jcm, agreements)
    parent.left.bp = new_jcm

    agreements = 0
    for i in range(len(parent.ancestral_sequence)):
        agreements += (p_seq[i] == '-' or r_seq[i] ==
                       '-' or p_seq[i] == r_seq[i])
    agreements /= len(parent.ancestral_sequence)
    new_jcm = np.full((4, 4), (1 - agreements) / 3)
    np.fill_diagonal(new_jcm, agreements)
    parent.right.bp = new_jcm


def train_model(orderings):
    for i in range(len(orderings)):
        for j in range(len(orderings[i])):
            if(orderings[i][j].left != None and orderings[i][j].right != None):
                set_ancestral_sequence(orderings[i][j])
                set_jcm_branch_lengths(orderings[i][j])


def train_initial_weightings(prob_matrix):
    probs = np.zeros(len(prob_matrix))
    total = 0
    for i in range(len(prob_matrix)):
        for j in range(len(prob_matrix[0])):
            probs[j] += prob_matrix[i][j]
            total += prob_matrix[i][j]
    return (probs / total)
