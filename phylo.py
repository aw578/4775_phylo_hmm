import math
import numpy as np
import nwk as parse_nwk


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


def build_orderings(init_models):
    final_orderings = []
    for i in range(len(init_models)):
        nwk = init_models[i]
        root, _ = parse_nwk.parse_nwk(nwk, 0)
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


def likelihood(ordering, data):
    base_conversion = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    m = len(data[0])

    for i in range(0, len(ordering)):
        current_node = ordering[i]
        name = current_node.name
        current_node.probs = np.zeros((m, 4))
        new_matrix = np.zeros((m, 4))

        # if it's a leaf node, 1s and 0s
        if(name is not None):
            for j in range(0, m):
                if(data[name][j] != '-'):
                    base_value = base_conversion[data[name][j]]
                    new_matrix[j][base_value] = 1
                else:
                    # should preserve values?
                    new_matrix[j][0] = 1
                    new_matrix[j][1] = 1
                    new_matrix[j][2] = 1
                    new_matrix[j][3] = 1

        # otherwise feisenstein
        else:
            l = current_node.left
            l_p = l.probs
            l_bp = l.bp
            r = current_node.right
            r_p = r.probs
            r_bp = r.bp

            for i in range(0, m):
                for x in range(0, 4):
                    l_sum = 0
                    r_sum = 0
                    for y in range(0, 4):
                        # might need to do log probs
                        l_sum += (l_bp[x][y] * l_p[i][y])
                        r_sum += (r_bp[x][y] * r_p[i][y])
                    new_matrix[i][x] = l_sum * r_sum
        # add back to probs
        current_node.probs = new_matrix
    # get current node matrix, add rows together, divide by 4, then add up logs (can reorder)
    final_matrix = ordering[len(ordering) - 1].probs
    likelihoods = np.zeros(m)
    sum = 0
    for i in range(0, m):
        for x in range(0, 4):
            likelihoods[i] += final_matrix[i][x]
        likelihoods[i] = np.log(likelihoods[i] / 4)
        sum += likelihoods[i]
    print(sum)
    return likelihoods


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    orderings: array of A orderings of trees
    data: array of N M-char sequences
Returns:
    likelihoods: A x M log likelihood matrix
'''


def phylo(orderings, data):
    a = len(orderings)
    m = len(data[0])
    likelihoods = np.zeros((a, m))
    for i in range(0, a):
        likelihoods[i] = likelihood(orderings[i], data)
    return likelihoods
