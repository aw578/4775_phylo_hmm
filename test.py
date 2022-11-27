import numpy as np


def read_models(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        trees = []
    for i in range(0, len(lines)):
        trees.append(lines[i].strip())
    return trees


read_models("hmm-sequence.fasta")
