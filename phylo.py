import numpy as np


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
