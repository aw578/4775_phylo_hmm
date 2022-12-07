import math
import numpy as np
import nwk


def jcm(len):
    at = 0.75*(1 - math.e**(-4 * len / 3))
    bp = np.full((4, 4), at / 3)
    for i in range(0, 4):
        bp[i][i] = (1 - at)
    return bp


def traverse(root, lst):
    if(root != None):
        traverse(root.left, lst)
        traverse(root.right, lst)
        lst.append(root)

# initialize ancestral sequences somewhere!


def build_orderings(init_models, obs):
    final_orderings = []
    for i in range(len(init_models)):
        model = init_models[i]
        root, _ = nwk.parse_nwk(model, 0, obs)
        ordering = []
        traverse(root, ordering)
        final_orderings.append(ordering)
    return final_orderings


def reweight(final_orderings, init_weights):
    for i in range(len(final_orderings)):
        for j in range(len(final_orderings[i])):
            final_orderings[i][j].branch_length *= init_weights[i]
            final_orderings[i][j].bp = jcm(final_orderings[i][j].branch_length)


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    data: sequence data (array of strings)
    models: trees (with branch lengths, names correspond to data pos)
    ordering = ordering of tree
Returns:
    likelihoods: 1 x M likelihood matrix
'''


def likelihood(ordering):
    m = len(ordering[0].ancestral_sequence)

    for i in range(0, len(ordering)):
        current_node = ordering[i]
        name = current_node.name

        # if it's a leaf node, 1s and 0s
        # TODO: can apparently be optimized
        if(name is not None):
            current_node.probs = np.zeros((m, 4))
            for j in range(0, m):
                if(ordering[i].ancestral_sequence[j] != 4):
                    current_node.probs[j][ordering[i].ancestral_sequence[j]] = 1
                else:
                    current_node.probs[j].fill(1)

        # otherwise feisenstein
        else:
            current_node.probs = np.multiply(
                np.dot(current_node.left.probs, current_node.left.bp), np.dot(current_node.right.probs, current_node.right.bp))

    # get current node matrix, add rows together, divide by 4, then add up logs (can reorder)
    final_matrix = ordering[len(ordering) - 1].probs
    likelihoods = np.log(final_matrix.sum(axis=1) / 4)
    # sum = np.sum(likelihoods)
    return likelihoods


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    orderings: array of A orderings of trees
    data: array of N M-char sequences
Returns:
    likelihoods: A x M log likelihood matrix
'''


def phylo(orderings):
    a = len(orderings)
    m = len(orderings[0][0].ancestral_sequence)
    likelihoods = np.zeros((a, m))
    for i in range(0, a):
        likelihoods[i] = likelihood(orderings[i])
    return likelihoods
