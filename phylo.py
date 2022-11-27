import math
import numpy as np


class Node():
    '''Calculates a branch probability matrix using the Jukes-Cantor model.
    '''

    def jcm(self, len):
        at = 0.75*(1 - math.e**(-4 * len / 3))
        bp = np.full((4, 4), at / 3)
        for i in range(0, 4):
            bp[i][i] = (1 - at)
        return bp
    ''' Initializes a node with given parameters.

    Arguments:
        name: name of node (corresponds to sequence number, so sequences should be left to right)
        left: left child (Node)
        right: right child (Node)
        branch_length: length of branch that leads to this node (float)
        probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None)
    '''

    def __init__(self, name, left, right, branch_length):
        self.name = name
        self.left = left
        self.right = right
        self.bp = self.jcm(branch_length)
        self.probs = None


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
    ordering = model.ordering

    for i in range(0, len(ordering)):
        current_node = ordering[i]
        name = current_node.name
        current_node.probs = np.zeros((m, 4))
        new_matrix = np.zeros((m, 4))

        # if it's a leaf node, 1s and 0s
        if(name is not None):
            for i in range(0, m):
                new_matrix[i][base_conversion[data[name][i]]] = 1

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
                        l_sum += (l_bp[x][y] * l_p[i][y])
                        r_sum += (r_bp[x][y] * r_p[i][y])
                    new_matrix[i][x] = l_sum * r_sum
        # add back to probs
        current_node.probs = new_matrix
    # get current node matrix, add rows together, divide by 4, then add up logs (can reorder)
    final_matrix = ordering[len(ordering) - 1].probs
    likelihoods = np.zeros(m)
    for i in range(0, m):
        for x in range(0, 4):
            likelihoods[i] += final_matrix[i][x]
        likelihoods[i] = np.log(likelihoods[i] / 4)

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
