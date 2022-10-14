# Sequence-Alingment
### Introduction
- this project is implemented to get sequence alignment[local, or global] for DNA, and Protein sequence. 

### Local Alignment
- It's also called Smith-Waterman Algorithm.

- If you have a long sequence and want to find any subsequences that are similar to any part of this sequence.

- Local alignment finds the best matching subsequences within two search sequences.

– Place a zero in any position in the scoring matrix if all of the other methods result in scores lower than zero.

### Global Alignment
- It's also called Needleman-Wunsch Algorithm

- N-W is guaranteed to find optimal alignments, although the algorithm does not search all possible alignments.

- It is an example of a dynamic programming algorithm:

- an optimal path (alignment) is identified by incrementallyextending optimal subpaths.

- Thus, a series of decisions is made at each step of the alignment to find the pair of residues with the best score. 

### Three steps to make alignment
1- set up a matrix.

2- score the matrix.

3- identify the optimal alignment(s).
