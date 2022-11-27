class Node():
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
        self.branch_length = branch_length
        self.bp = None
        self.probs = None


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


def parse_leaf(nwk, index):
    end_index = index
    while(nwk[end_index].isnumeric()):
        end_index += 1
    name = int(nwk[index:end_index])
    index = end_index

    # get past colon
    index += 1
    branch_length, end_index = parse_branch_length(nwk, index)
    return Node(name, None, None, branch_length), end_index


def parse_nwk(nwk, index):
    index += 1
    if(nwk[index].isnumeric()):
        left, index = parse_leaf(nwk, index)
        index += 2
        right, index = parse_leaf(nwk, index)
        index += 2
        branch_length, end_index = parse_branch_length(nwk, index)
        return Node(None, left, right, branch_length), end_index
    else:
        left, index = parse_nwk(nwk, index)
        index += 2
        right, index = parse_nwk(nwk, index)
        index += 2
        branch_length, end_index = parse_branch_length(nwk, index)
        return Node(None, left, right, branch_length), end_index
