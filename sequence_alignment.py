### Importing Packages ###
from Bio.SubsMat import MatrixInfo
import Blosum62Matrix
import numpy as np
import pandas as pd
b62 = MatrixInfo.blosum62


# this function take which sequence the user want to make alignment for it
def SequenceAlignment(selectedSequence):
    # DNA Sequence
    if selectedSequence == 1:
        UserSelect = int(input(
            "***** DNA Sequence *****\n Enter The Number Of Your Choice.\n  1.Local Alignment \n  2.Global Alignment\n"))
        if UserSelect == 1:
            DNA_LocalAlignment()
        elif UserSelect == 2:
            DNA_GlobalAlignment()

    # Protein Sequence
    elif selectedSequence == 2:
        UserSelect = int(input(
            "***** Protein Sequence *****\n Enter The Number Of Your Choice.\n  1.Local Alignment \n  2.Global Alignment\n"))
        if UserSelect == 1:
            Protein_LocalAlignment()
        elif UserSelect == 2:
            Protein_GlobalAlignment()

    else:
        print("Please Enter a Valid Choice.")


# this function take which DNA sequence the user want to make local alignment for it
def DNA_LocalAlignment():
    Seqence1 = str.upper(input("Enter The First DNA Sequence: "))
    Seqence2 = str.upper(input("Enter The Second DNA Sequence: "))

    DNAseq1Length = len(Seqence1)
    DNAseq2Length = len(Seqence2)

    DNALoopMatrix = [[0 for x in range(DNAseq1Length + 1)] for y in range(DNAseq2Length + 1)]
    for i in range(DNAseq2Length + 1):
        DNALoopMatrix[i][0] = 0
    for j in range(DNAseq1Length + 1):
        DNALoopMatrix[0][j] = 0
    for i in range(1, DNAseq2Length + 1, 1):
        for j in range(1, DNAseq1Length + 1, 1):
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = 1
            else:
                MatrixScore = -2
            DNALoopMatrix[i][j] = max(DNALoopMatrix[i - 1][j] - 1, DNALoopMatrix[i][j - 1] - 1,
                                      DNALoopMatrix[i - 1][j - 1] + MatrixScore, 0)
    #print(DNALoopMatrix[DNAseq2Length][DNAseq1Length])

    Seq1Result = ""
    Seq2Result = ""
    SeqMatch = ""
    i = DNAseq2Length
    j = DNAseq1Length
    while i > 0 and j > 0:
        up = DNALoopMatrix[i - 1][j] - 1
        left = DNALoopMatrix[i][j - 1] - 1
        if Seqence1[j - 1] == Seqence2[i - 1]:
            MatrixScore = 1
        else:
            MatrixScore = -2
        MatrixDiagonal = DNALoopMatrix[i - 1][j - 1] + MatrixScore
        if DNALoopMatrix[i][j] == MatrixDiagonal:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += Seqence2[i - 1]
            if MatrixScore == 1:
                SeqMatch += "|"
            else:
                SeqMatch += " "
            i -= 1
            j -= 1
        elif DNALoopMatrix[i][j] == up:
            Seq1Result += "-"
            SeqMatch += " "
            Seq2Result += Seqence2[i - 1]
            i -= 1
        else:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

    while (i > 0):
        Seq1Result += "-"
        Seq2Result += Seqence2[i - 1]
        SeqMatch += " "
        i -= 1
    while (j > 0):
        Seq1Result += Seqence1[j - 1]
        Seq2Result += "-"
        SeqMatch += " "
        j -= 1

    Seq1Result = Seq1Result[::-1]
    SeqMatch = SeqMatch[::-1]
    Seq2Result = Seq2Result[::-1]

    print("\n*** DNA Local Alignment Output ***")
    print(Seq1Result)
    print(SeqMatch)
    print(Seq2Result)


# this function take which DNA sequence the user want to make global alignment for it
def DNA_GlobalAlignment():
    Seqence1 = str.upper(input("Enter The First DNA Sequence: "))
    Seqence2 = str.upper(input("Enter The Second DNA Sequence: "))

    DNAseq1Length = len(Seqence1)
    DNAseq2Length = len(Seqence2)
    DNALoopMatrix = [[0 for x in range(DNAseq1Length + 1)] for y in range(DNAseq2Length + 1)]
    for i in range(DNAseq2Length + 1):
        DNALoopMatrix[i][0] = -i
    for j in range(DNAseq1Length + 1):
        DNALoopMatrix[0][j] = -j
    for i in range(1, DNAseq2Length + 1, 1):
        for j in range(1, DNAseq1Length + 1, 1):
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = 1
            else:
                MatrixScore = -2
            DNALoopMatrix[i][j] = max(DNALoopMatrix[i - 1][j] - 1, DNALoopMatrix[i][j - 1] - 1,
                                      DNALoopMatrix[i - 1][j - 1] + MatrixScore)
    #print(DNALoopMatrix[DNAseq2Length][DNAseq1Length])

    Seq1Result = ""
    Seq2Result = ""
    SeqMatch = ""
    i = DNAseq2Length
    j = DNAseq1Length
    while i > 0 and j > 0:
        up = DNALoopMatrix[i - 1][j] - 1
        left = DNALoopMatrix[i][j - 1] - 1
        if Seqence1[j - 1] == Seqence2[i - 1]:
            MatrixScore = 1
        else:
            MatrixScore = -2
        diagonal = DNALoopMatrix[i - 1][j - 1] + MatrixScore
        if DNALoopMatrix[i][j] == diagonal:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += Seqence2[i - 1]
            if MatrixScore == 1:
                SeqMatch += "|"
            else:
                SeqMatch += " "
            i -= 1
            j -= 1
        elif DNALoopMatrix[i][j] == up:
            Seq1Result += "-"
            SeqMatch += " "
            Seq2Result += Seqence2[i - 1]
            i -= 1
        else:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

    while (i > 0):
        Seq1Result += "-"
        Seq2Result += Seqence2[i - 1]
        SeqMatch += " "
        i -= 1
    while (j > 0):
        Seq1Result += Seqence1[j - 1]
        Seq2Result += "-"
        SeqMatch += " "
        j -= 1

    Seq1Result = Seq1Result[::-1]
    SeqMatch = SeqMatch[::-1]
    Seq2Result = Seq2Result[::-1]

    print("\n*** DNA Global Alignment Output ***")

    print(Seq1Result)
    print(SeqMatch)
    print(Seq2Result)



# this function take which Protein sequence the user want to make global alignment for it
def Protein_GlobalAlignment():
    # Sequence1 = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFK'
    # Sequence2 = 'MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKGHGKKVVAALIEAANHIDDIAGTLSKLSDLHAHKLRVDPVNFK'
    Seqence1 = str.upper(input("Enter The First Protein Sequence: "))
    Seqence2 = str.upper(input("Enter The Second Protein Sequence: "))

    ProteinGapScore = -1
    Proteinseq1Length = len(Seqence1)
    Proteinseq2Length = len(Seqence2)
    ProteinLoopMatrix = [[0 for x in range(Proteinseq1Length + 1)] for y in range(Proteinseq2Length + 1)]
    for i in range(Proteinseq2Length + 1):
        ProteinLoopMatrix[i][0] = -i
    for j in range(Proteinseq1Length + 1):
        ProteinLoopMatrix[0][j] = -j
    for i in range(1, Proteinseq2Length + 1, 1):
        for j in range(1, Proteinseq1Length + 1, 1):
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
                # print('M= ', Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]))
            else:
                MatrixScore = Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
                # print('MM= ', Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]))
            ProteinLoopMatrix[i][j] = max(ProteinLoopMatrix[i - 1][j] + ProteinGapScore,
                                          ProteinLoopMatrix[i][j - 1] + ProteinGapScore,
                                          ProteinLoopMatrix[i - 1][j - 1] + int(MatrixScore))
    print(ProteinLoopMatrix[Proteinseq2Length][Proteinseq1Length])

    Seq1Result = ""
    Seq2Result = ""
    SeqMatch = ""
    i = Proteinseq2Length
    j = Proteinseq1Length

    while i > 0 and j > 0:
        up = ProteinLoopMatrix[i - 1][j] - 1
        left = ProteinLoopMatrix[i][j - 1] - 1
        if Seqence1[j - 1] == Seqence2[i - 1]:
            MatrixScore = Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
            # print('M= ', Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]))
        else:
            MatrixScore = Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
            # print('MM= ', Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]))

        diagonal = ProteinLoopMatrix[i - 1][j - 1] + MatrixScore
        if ProteinLoopMatrix[i][j] == diagonal:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += Seqence2[i - 1]
            if MatrixScore == Blosum62Matrix.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]):
                SeqMatch += "|"
            else:
                SeqMatch += " "
            i -= 1
            j -= 1
        elif ProteinLoopMatrix[i][j] == up + ProteinGapScore:
            Seq1Result += "-"
            SeqMatch += " "
            Seq2Result += Seqence2[i - 1]
            i -= 1
        else:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

    while (i > 0):
        Seq1Result += "-"
        Seq2Result += Seqence2[i - 1]
        SeqMatch += " "
        i -= 1
    while (j > 0):
        Seq1Result += Seqence1[j - 1]
        Seq2Result += "-"
        SeqMatch += " "
        j -= 1

    Seq1Result = Seq1Result[::-1]
    SeqMatch = SeqMatch[::-1]
    Seq2Result = Seq2Result[::-1]

    print("\n*** Protein Global Alignment Output ***")
    print(Seq1Result)
    print(SeqMatch)
    print(Seq2Result)


# this function take which Protein sequence the user want to make global alignment for it
def Protein_LocalAlignment():
    Seqence1 = str.upper(input("Enter The First Protein Sequence: "))
    Seqence2 = str.upper(input("Enter The Second Protein Sequence: "))

    Proteinseq1Length = len(Seqence1)
    Proteinseq2Length = len(Seqence2)
    ProteinLoopMatrix = [[0 for x in range(Proteinseq1Length + 1)] for y in range(Proteinseq2Length + 1)]

    i = 0
    j = 0
    ProteinGapScore = -1
    for i in range(1, Proteinseq2Length + 1, 1):
        for j in range(1, Proteinseq1Length + 1, 1):
            if Seqence1[j - 1] == Seqence2[i - 1]:
                MatrixScore = Blosum62.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
            else:
                MatrixScore = Blosum62.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
            ProteinLoopMatrix[i][j] = max(ProteinLoopMatrix[i - 1][j] + ProteinGapScore,
                                          ProteinLoopMatrix[i][j - 1] + ProteinGapScore,
                                          ProteinLoopMatrix[i - 1][j - 1] + MatrixScore, 0)
    #print(ProteinLoopMatrix[Proteinseq2Length][Proteinseq1Length])

    Seq1Result = ""
    Seq2Result = ""
    SeqMatch = ""
    i = Proteinseq2Length
    j = Proteinseq1Length
    while i > 0 and j > 0:
        up = ProteinLoopMatrix[i - 1][j] - 1
        left = ProteinLoopMatrix[i][j - 1] - 1
        if Seqence1[j - 1] == Seqence2[i - 1]:
            MatrixScore = Blosum62.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
        else:
            MatrixScore = Blosum62.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1])
        MatrixDiagonal = ProteinLoopMatrix[i - 1][j - 1] + MatrixScore
        if ProteinLoopMatrix[i][j] == MatrixDiagonal:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += Seqence2[i - 1]
            if MatrixScore == Blosum62.substitution_matrix(Seqence1[j - 1], Seqence2[i - 1]):
                SeqMatch += "|"
            else:
                SeqMatch += " "
            i -= 1
            j -= 1
        elif ProteinLoopMatrix[i][j] == up + ProteinGapScore:
            Seq1Result += "-"
            SeqMatch += " "
            Seq2Result += Seqence2[i - 1]
            i -= 1
        else:
            Seq1Result += Seqence1[j - 1]
            Seq2Result += "-"
            SeqMatch += " "
            j -= 1

    while (i > 0):
        Seq1Result += "-"
        Seq2Result += Seqence2[i - 1]
        SeqMatch += " "
        i -= 1
    while (j > 0):
        Seq1Result += Seqence1[j - 1]
        Seq2Result += "-"
        SeqMatch += " "
        j -= 1

    Seq1Result = Seq1Result[::-1]
    SeqMatch = SeqMatch[::-1]
    Seq2Result = Seq2Result[::-1]

    print("\n*** Protein Local Alignment Output ***")
    print(Seq1Result)
    print(SeqMatch)
    print(Seq2Result)


### Main ###
selectedSequence = int(input("***** Sequence Alignment *****\n Enter The Number Of Your Choice.\n  1.DNA Sequence \n  2.Protein Sequence\n"))
SequenceAlignment(selectedSequence)
