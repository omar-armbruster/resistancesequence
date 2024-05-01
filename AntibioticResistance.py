import numpy as np
import sys
from Bio.Align import substitution_matrices

def create_scoring_matrix(m, n, sigma):
    # Make all nodes 0 b/c local alignment can start from anywhere
    return np.zeros((m+1, n+1), dtype=int)


def fill_scoring_matrix(x, y, S, sigma, scoring_matrix):
    # Fill scoring matrix
    m, n = len(x), len(y)
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            #scores for possible alignments
            insert = S[i, j-1] - sigma
            delete = S[i-1, j] - sigma
            match = S[i-1, j-1] + scoring_matrix[(x[i-1], y[j-1])]

            S[i, j] = max(match, delete, insert, 0) #recurrence relation

            #keep track
            if S[i, j] > max_score:
                max_score = S[i, j]
                max_pos = (i, j)

    return S, max_score, max_pos


def backtrack(S, x, y, start_pos, sigma, scoring_matrix):
    #Backtrack to find optimal local alignment
    X_align, Y_align = "", ""  # Initialize aligned strings
    i, j = start_pos  # Start backtracking from the position of the max score

    # Backtrack until reaching a cell with score 0 (sigals end of local alignment)
    while i > 0 and j > 0 and S[i, j] > 0:
        if S[i, j] == S[i-1, j-1] + scoring_matrix[(x[i-1], y[j-1])]:
            # Diagonal: match/mismatch
            X_align = x[i-1] + X_align
            Y_align = y[j-1] + Y_align
            i -= 1
            j -= 1
        elif S[i, j] == S[i, j-1] - sigma:
            #Left: insertion
            X_align = "-" + X_align
            Y_align = y[j-1] + Y_align
            j -= 1
        else: # S[i, j] == S[i-1, j] - sigma:
            #Up: deletion
            X_align = x[i-1] + X_align
            Y_align = "-" + Y_align
            i -= 1

    return X_align, Y_align


def local_alignment(x, y, sigma=5):
    #Perform local alignment using PAM250 scoring matrix
    scoring_matrix = substitution_matrices.load("PAM250")
    m, n = len(x), len(y)

    S = create_scoring_matrix(m, n, sigma)
    S, max_score, max_pos = fill_scoring_matrix(x, y, S, sigma, scoring_matrix)

    X_align, Y_align = backtrack(S, x, y, max_pos, sigma, scoring_matrix)


    return max_score, X_align, Y_align


# def parse_file(filename):
#     sequences = []
#     with open(filename, 'r') as file:
#         for line in file:
#             sequence = line.strip()
#             sequences.append(sequence)
#     return sequences


def pairwiseAlign(genomes):
    #Aligns each genome with all of the other genomes. Puts them in list rank and then sorts
    #based on the max score. List is returned with elements of [score, xalign, yalign]
    rank = []
    for i in range((len(genomes))):
        for j in range(i+1, len(genomes)):
          align = local_alignment(genomes[i], genomes[j], 5)
          rank.append(align)
    rank = sorted(rank, key=lambda x : x[0], reverse = True)
    return rank

def createProfile(alignment):
    x = alignment[1]
    y = alignment[2]
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    prof = np.zeros((len(amino_acids), len(x)))
    for i in range(len(x)):
        if x[i] == y[i]:
            prof[amino_acids.index(x[i])][i] += 1
    
    return prof
    
    

if __name__ == "__main__":
    sequences = ["MKTIIALSYIFCLV","TIIALSYIFCLVFA","ALSYIFCLVFADYK","CLVFADYKDDDDK","IFCLVFADY","SYIFCLVFA"]
    all_alignments = pairwiseAlign(sequences)
    two_seq = all_alignments[0]
    profile = createProfile(two_seq)
    print(profile)


print(substitution_matrices.load("PAM250"))


# if __name__ == "__main__":
#     filename = sys.argv[1]
#     sequences = parse_file(filename)
#     x = sequences[0]
#     y = sequences[1]
#     max_score, X_align, Y_align = local_alignment(x, y)
#     print(max_score)
#     print(X_align)
#     print(Y_align)