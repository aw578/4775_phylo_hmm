import numpy as np


class Node():
    ''' Initializes a node with given parameters.

    Arguments:
        name: name of node (corresponds to sequence number, so sequences should be left to right)
        left: left child (Node)
        right: right child (Node)
        branch_length: length of branch that leads to this node (float)
        probs: probability of observed bases beneath this node
                [4 x m probs for 'ACGT'] (initialized to None)
    '''

    def __init__(self, name, left, right, branch_length, ancestral_sequence):
        self.name = name
        self.left = left
        self.right = right
        self.branch_length = branch_length
        self.ancestral_sequence = ancestral_sequence
        self.bp = None
        self.probs = np.zeros(len(ancestral_sequence))


def parse_branch_length(nwk, index):
    end_index = index
    while(end_index < len(nwk) and (nwk[end_index].isnumeric() or nwk[end_index] == '.' or nwk[end_index] == "-")):
        end_index += 1
    if(end_index == index):
        # root
        branch_length = 0
    else:
        branch_length = float(nwk[index:end_index])
    return branch_length, end_index


def parse_leaf(nwk, index, obs):
    end_index = index
    while(nwk[end_index].isnumeric()):
        end_index += 1
    name = int(nwk[index:end_index])
    index = end_index

    # get past colon
    index += 1
    branch_length, end_index = parse_branch_length(nwk, index)
    return Node(name, None, None, branch_length, obs[name]), end_index


def parse_nwk(nwk, index, obs):
    index += 1
    if(nwk[index].isnumeric()):
        left, index = parse_leaf(nwk, index, obs)
        index += 2
        if(nwk[index].isnumeric()):
            right, index = parse_leaf(nwk, index, obs)
        else:
            right, index = parse_nwk(nwk, index, obs)
        index += 2
        branch_length, end_index = parse_branch_length(nwk, index)
        return Node(None, left, right, branch_length, np.zeros(len(obs[0]))), end_index
    else:
        left, index = parse_nwk(nwk, index, obs)
        index += 2
        if(nwk[index].isnumeric()):
            right, index = parse_leaf(nwk, index, obs)
        else:
            right, index = parse_nwk(nwk, index, obs)
        index += 2
        branch_length, end_index = parse_branch_length(nwk, index)
        return Node(None, left, right, branch_length, np.zeros(len(obs[0]))), end_index
