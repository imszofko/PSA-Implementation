import requests, sys
from Bio import SeqIO

#Establishing the API to ENSMBL
#rest.ensembl.org/sequence/id/ENSG00000139618                                 ##This is the direct link to the sequence

#Keeping both inputs for data files together so they are entered one after the other
ensemblID = input(str("Enter ENSEMBL Sequence ID: "))
inputFile = input(str("Enter the pathway and file name of the saved fasta data: "))

server = "https://rest.ensembl.org"
extension = f"/sequence/id/{ensemblID}"

ensemblResponse = requests.get(server + extension, headers = {"Content-Type" : "text/x-fasta"})                 #Fetched sequence stored here

if ensemblResponse.ok:
    # Get FASTA content as a string
    fasta_content = ensemblResponse.text
    print("Fasta content saved")
    # Extract the sequence (ignoring the header)
    fasta_lines = fasta_content.strip().splitlines()
    header = fasta_lines[0]  # First line is the header
    sequence = ''.join(fasta_lines[1:])  # Remaining lines are the sequence
    print("Sequence Extracted.")

    # Save to a file
    with open("Seq_ENSEMBL.fasta", "w") as file:
        file.write(fasta_content)  # Save as is
    print("Sequence written to file.")
if not ensemblResponse.ok:
  ensemblResponse.raise_for_status()
  sys.exit()

#Global variable for FASTA file from ENSEMBL
globalFastaFile = "Seq_ENSEMBL.fasta"                               #Integrated into the parameter since this file name isnt changing

def fastaFile(globalFasta, inputFastaFile):
    #read the global FASTA data file
    globalDataFileseq1 = list(SeqIO.parse(globalFasta, "fasta"))
    
    #read the input FASTA data file
    variantDataseq2 = list(SeqIO.parse(inputFastaFile, "fasta"))
    
    #Extracting IDs and sequences for the first sequence
    for record1 in globalDataFileseq1:
        seq1 = record1.seq
        seq1ID = record1.id
        print("Sequence one:" , seq1ID + '\n' + seq1[0:10])
    #Extracting IDs and sequences for the second sequence
    for record2 in variantDataseq2:
        seq2 = record2.seq
        seq2ID = record2.id
        print("Sequence two:", seq2ID + '\n' + seq2[0:10])

    return seq1, seq1ID, seq2, seq2ID
seq1, seq1ID, seq2, seq2ID = fastaFile(globalFastaFile, inputFile)

"""
"Testing with these variables instead of big data files
seq1 = str("actagactctcagacgtacgatcgatcgagcatctagcatcgatcCTAGACGCTATCTAGCGAGCACTTACGGCGGCTCTCGCGCGCGGCGCATTATATCTCTGAATCTGCGTAATCGGCTAGCTAGCTAGCAGAGCTagcatgctacgatcatcgatcatcgatcgatcatcgactaa").upper()
seq1ID = "Name of Seq1"
seq2 = str("actgacgactgatcgagcatctagcatcgacgactatcgatcgagcatctagcatcggactatctcttctactattatgatcgagcatctagcatcgtttgatcgagcatctagcatcgtttatcgagcatcgatcagagtcgatcattcagctcagATCGGCTAGCTTCGATTTCTaa").upper()
seq2ID = "Name of Seq1"
"""

def smithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0):
    #Create the matrix
    s1 = len(seq1)                  #Sequence one is the vertical sequence to the left
    s2 = len(seq2)                  #Sequence two is teh horizontal sequence on top
    
    tempMat = []
    for i in range(s1 + 1):
        tempMat.append([])
        for j in range(s2 + 1):
            tempMat[-1].append(0)   #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
    print("Matrix made.")

    print("Starting to fill matrix.")
    for i in range(1, s1 + 1):
        for j in range(1, s2 + 1):
            if seq1[i - 1] == seq2[j - 1]:                              ##If the diagonals are equal then it is a match
                scoreMatch = tempMat[i - 1][j - 1] + match                  
            else:
                scoreMatch = tempMat[i - 1][j - 1] + mismatch           ##If not then it is a mismatch
            
            scoreInsert = tempMat[i][j - 1] + gap                       ##How the cell is filled in if it is a gap to the left
            scoreDelete = tempMat[i - 1][j] + gap                       ##Gap down

            tempMat[i][j] = max(0, scoreMatch, scoreInsert, scoreDelete)       ##Searching for the max score from the three possibilities

            if tempMat[i][j] >= maxScore:
                maxScore = tempMat[i][j]
                mismatch, gap = (i, j)
    print("Matrix filled.")

    print("Starting traceback.")
    #Traceback starts at the element with the highest score
    #Finding the element with the highest score
    alignA = ""
    alignB = ""
    maximum = 0
    for row in range(1, s1 + 1):
        for column in range(1, s2 + 1):
            if maximum < tempMat[row][column]:
                maximum = tempMat[row][column]           #max value is stored
                i = row                                     #Storing the row 
                j = column                                  #Storing column

    #Backtracing
    '''
    Based on the source of each score recursively until 0 is encountered
    Segments that have the highest similarity score based on the given scoring system is generated in this process
    to obtain the second best local alignment, apply the traceback process starting at the second highest score outside the trace of the best alignment
    Very important is that no negative score is assigned in the scoring system which enables local alignment
    '''
    while tempMat[i][j] != 0:
        if seq1[i - 1] == seq2[j - 1]:
            alignA += seq1[i - 1]
            alignB += seq2[j - 1]
            i -= 1
            j -= 1
        #This is for if the sequences dont match
        elif seq1[i - 1] != seq2[j - 1]:
            tempList = [tempMat[i - 1][j - 1], tempMat[i - 1][j], tempMat[i][j - 1]]
            if max(tempList) == tempList[0]:
                alignA += seq1[i - 1]
                alignB += seq2[j - 1]
                i -= 1
                j -= 1
            #If the maximum value is the 1st indexed position, topvalue
            elif max(tempList) == tempList[1]:
                alignA += seq1[i - 1]
                alignB += '-'
                i -= 1
            #If the max value is the second indexed position, leftvalue
            elif max(tempList) == tempList[-1]:
                alignA += '-'
                alignB += seq2[j - 1]
                j -= 1
        else:
            print('Error. Exit.')
            i = 0
            j = 0

    print("Backtrace complete")

    #Reverse the strings
    alignA = alignA[::-1]
    alignB = alignB[::-1]
    print("Sequences reversed")

    #Storing match  and such  scores and symbols
    matchString = ""
    for i in range (len(alignA)):
        if alignA[i] == alignB[i]:          #When A = B then line is appended to matchString
            matchString += '|'
        elif alignA[i] != alignB[i]:
            if (alignA[i] == '-' or alignB == '-'):
                matchString += " "          #When not a match then blank is appended
            else:
                matchString += "*"          #Or an asterisk
    print("Match scores noted.")

    ##Calculating alignment scores of the local alignment
    alignScore = 0
    for i in range(len(matchString)):       ##Calculated based on the matchString.
        if matchString[i] == '|':           ##Lines are a +1 and asterisk or empty space is -1
            alignScore += 1                 ##Assigned value is then appended to alignScore
        elif (matchString[i] == '*' or matchString[i] == " "):
            alignScore += -1
    print("Alignment score calculated.")
    return alignA, alignB, matchString, alignScore

#Printing out the results
output1, output2, matchString, alignScore = smithWaterman(seq1, seq2, mismatch = -1, gap = -1, match = 1, maxScore = 0)
print("Alignment Score: ", alignScore)

#To write the sequences to a .txt file to look like EMBOSS and 'highlight' nucleotides that are different
'''
##Adding the highlighting function. When there are two or more consercutive nucleotides that are not 
##matching between the two sequences, they will be represeneted in lowercase letters in the txt file.
##The txt file should be fluid and easy to read. It will have the seq1 and seq2 and its match string broken into chunks for easy reading.
##The difference will be in lowercase. And the sequence name will be stated in the beginnning of each sequence, even in the chunks
'''
def highlightRegions(output1, output2):
    print("Starting to analyze and highlight regions with differences.")
    highlightSeq1 = []
    highlightSeq2 = []

    for out1, out2 in zip(output1, output2):
        if out1 != out2 : #mismatch or gap
            # Mismatches are converted to lowercase
            highlightSeq1.append(out1.lower() if out1 != '-' else '-')
            highlightSeq2.append(out2.lower() if out2 != '-' else '-')
        else: #match
            # Matches and gaps are kept as-is
            highlightSeq1.append(out1)
            highlightSeq2.append(out2)

    return ''.join(highlightSeq1), ''.join(highlightSeq2) ##Joins the strings together
print("Finished highlighting.")

def writeToFile(output1, output2, matchString, seq1ID, seq2ID, output_file):
    print("Started writing alignment to text file.")
    #calls up highlightRegion function and stores the alignments in 
    highlightSeq1, highlightSeq2 = highlightRegions(output1, output2)
    chunk_size = 60                                                     ##How wide the chunks of the alignment will be displayed as
    ##Had big issue with the match string going in between the sequence IDs, this is a good fix
    seqPadding = max(len(seq1ID), len(seq2ID)) + 2                      ##To make sure that the matchstring starts and is accurate to the alignment
    with open(output_file, "w") as file:
        for start in range(0, len(output1), chunk_size):                ##Formatting how everything should be written into the file
            chunk1 = highlightSeq1[start:start + chunk_size]
            chunk2 = highlightSeq2[start:start + chunk_size]
            matchChunk = matchString[start:start + chunk_size]

            file.write(f"{seq1ID.ljust(seqPadding)}{chunk1}\n")
            file.write(f"{' ' * seqPadding}{matchChunk}\n")
            file.write(f"{seq2ID.ljust(seqPadding)}{chunk2}\n\n")

    print(f"Highlighted alignment written to {output_file}")

writeToFile(output1, output2, matchString, seq1ID, seq2ID, output_file = "SW_Alignment.txt")
##End.