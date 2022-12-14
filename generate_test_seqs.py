import random

import numpy as np
import phylo

seq = [random.randrange(0, 4, 1) for _ in range(10000)]

conserved = [phylo.jcm(0.001), phylo.jcm(0.0002)]
divergent = [phylo.jcm(0.1), phylo.jcm(0.002)]

final_seqs = np.zeros((2, 100000), int)

for i in range(0, 2000):
    final_seqs[0][i] = np.random.choice(4, 1, p=conserved[0][seq[i]])
    final_seqs[1][i] = np.random.choice(4, 1, p=conserved[1][seq[i]])
for i in range(2000, 4000):
    final_seqs[0][i] = np.random.choice(4, 1, p=divergent[0][seq[i]])
    final_seqs[1][i] = np.random.choice(4, 1, p=divergent[1][seq[i]])
for i in range(4000, 6000):
    final_seqs[0][i] = np.random.choice(4, 1, p=conserved[0][seq[i]])
    final_seqs[1][i] = np.random.choice(4, 1, p=conserved[1][seq[i]])
for i in range(6000, 8000):
    final_seqs[0][i] = np.random.choice(4, 1, p=divergent[0][seq[i]])
    final_seqs[1][i] = np.random.choice(4, 1, p=divergent[1][seq[i]])
for i in range(8000, 10000):
    final_seqs[0][i] = np.random.choice(4, 1, p=conserved[0][seq[i]])
    final_seqs[1][i] = np.random.choice(4, 1, p=conserved[1][seq[i]])
with open("test_bcdn.fasta", "w") as f:
    i = 0
    while(i < 10000):
        j = 0
        while(j < 60 and i < 10000):
            if(i >= 0 and i < 2000):
                f.write("C")
            elif(i >= 2000 and i < 4000):
                f.write("D")
            elif(i >= 4000 and i < 6000):
                f.write("C")
            elif(i >= 6000 and i < 8000):
                f.write("D")
            else:
                f.write("C")
            j += 1
            i += 1
        f.write('\n')

with open("test_seq.fasta", "w") as f:
    f.write("> sequence 0")
    f.write("\n")
    i = 0
    while(i < 10000):
        j = 0
        while(j < 60 and i < 10000):
            if(final_seqs[0][i] == 0):
                f.write('A')
            elif(final_seqs[0][i] == 1):
                f.write('C')
            elif(final_seqs[0][i] == 2):
                f.write('G')
            else:
                f.write('T')
            j += 1
            i += 1
        f.write('\n')

    f.write("> sequence 1")
    f.write("\n")
    i = 0
    while(i < 10000):
        j = 0
        while(j < 60 and i < 10000):
            if(final_seqs[1][i] == 0):
                f.write('A')
            elif(final_seqs[1][i] == 1):
                f.write('C')
            elif(final_seqs[1][i] == 2):
                f.write('G')
            else:
                f.write('T')
            j += 1
            i += 1
        f.write('\n')
