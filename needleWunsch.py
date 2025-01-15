import requests, sys
from Bio import SeqIO


#Establishing the API to ENSMBL
##rest.ensembl.org/sequence/id/ENSG00000139618                     #Link to put into web browser that will bring you to the sequence for BRCA2 gene

#Keeping both inputs for data files together so they are entered one after the other
ensemblID = input(str("Enter ENSEMBL Sequence ID: "))
inputFile = input(str("Enter the pathway and file name of the saved fasta data: "))

server = "https://rest.ensembl.org"
extension = f"/sequence/id/{ensemblID}"

ensemblResponse = requests.get(server + extension, headers = {"Content-Type" : "text/x-fasta"})             #Fetched sequence stored here

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

#Defining global variable for FASTA file path
#Global fasta file is always the same since the file name is created as the saved content of the fetched data 
globalFastaFile = "Seq_ENSEMBL.fasta"

def FastaFile(globalFasta, inputFastaFile):
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
seq1, seq1ID, seq2, seq2ID = FastaFile(globalFastaFile, inputFile)


'''
#Testing with these variables instead of big data files
seq1 = str("ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCACAACGTCAacGACAACGTCAacGACAACGTCAacGACAACAACGTCAacGACGTCAacGACAACGTCAacGACAACACAACGTCAacGGTCAacGACAACGTCAacGACAACGTCAacGCTATGG").upper()
seq1ID = "Name of Seq1"
seq2 = str("CTGTCTCTGCGACACAACGTCAacGAGACAACGTCAacGacaACACAACGTCAacGACAacGCGTCAacGgatT").upper()
seq2ID = "Name of Seq1"
'''

##Matrix of Zeros Function
def Zeros(rows, cols):
    #Define an empty list
    matZero = []
    #Set up rows in the matrix
    for x in range(rows):
        #For each row add an empty list
        #this appends an empty array for each x in rows
        matZero.append([])
        #Set up the columns in each row
        for y in range(cols):
            #add zero to each y in cols in each row
            #-1 index is always the last index in an array, it will add 0s to the the last index of the original adding arrays and also the 0s in.
            matZero[-1].append(0)
    return matZero

#Return the score between any two bases in alignment
def MatchScore(match1, match2, mismatch = -1, gap = -1, match = 1):
    if match1 == match2:
        return match                #+1
    elif match1 == '-' or match2 == '-':
        return gap                  #-1
    else:
        return mismatch             #-1
    
##Function that fills out the matrix of scores
def NeedleWunsch(seq1, seq2, gap = -1, match = 1):
    print("Starting NW sequence alignment.")
    
    #Length of two sequence
    s1 = len(seq1) #n
    s2 = len(seq2) #m

    #FILLING IN SCOREMATRIX

    #Generate the matrix of Zeros to stores the scores
    scoreMatrix = Zeros(s1 + 1, s2 + 1)                     #+1 to create the 0s col and rows
    print("Starting to make scoring matrix.")

    #Calculate score table
    #First col and first row of Zeros
    
    #First column
    for i in range(0, s1 + 1):
        scoreMatrix[i][0] = gap * i

    #Rows
    for j in range(0, s2 + 1):
        scoreMatrix[0][j] = gap * j

    for i in range(1, s1 + 1):
        for j in range(1, s2 + 1):

            #Calculating the score by checking the top, the left, and the diagonal squares
            
            match = scoreMatrix[i - 1][j - 1] + MatchScore(seq1[i - 1], seq2[j - 1])                  #Diagonal one from [i][j] position at that time
            insertGap = scoreMatrix[i - 1][j] + gap                                               #Up one of [i][j] position at that time
            misDelete = scoreMatrix[i][j - 1] + gap                                               #To the left of [i][j] position at that time

            #Now going to recording the maximum score from the three possibilities calculated
            scoreMatrix[i][j] = max(match, insertGap, misDelete)
    print("Finished filling matrix.")

    print("Starting traceback.")
    #Traceback from Wikipedia pseudocode
    #Creating a variable to store each alignment
    alignA = ""
    alignB = ""

    ##Storing the length of both sequences
    i = s1 
    j = s2

    ##While loop to trace where we are in the matrix during the traceback
    ##The traceback starts in the lower right square of the matrix and moves backwards throughout the matrix
    ##Starts in bottom right corner and compares the value with the three possible sources to see which it comes from
    #If seq1 and seq2 are aligned they match and it is +1 and the same with the rest of the stuff essentially

    while i > 0 and j > 0:
        #checking that it is a match, if it is then it appends to alignA/B and then jump to the diagnol value after
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + MatchScore(seq1[i - 1], seq2[j - 1]):                           #Diagonal step
            alignA += seq1[i - 1]
            alignB += seq2[j - 1]
            i -= 1 #move forward
            j -= 1
        
        #Checking of the index [i][j] and one above are equal with and to update i and j to correspond to that cell
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap:                                                          #Step above the [i][j]
            alignA += seq1[i - 1]
            alignB += '-'
            i -= 1
        
        elif scoreMatrix[i][j] == scoreMatrix[i][j - 1] + gap:                                                          #Step to the left of [i][j]
            alignA += '-'
            alignB += seq2[j - 1]
            j -= 1

    #WIll append to alignA or B depending on the result of j or i (seq 1 or seq2) and finish tracing up to the top left cell
    ##Appending occurs also before this line (ln 106 to 123)
    
    while i > 0:
        alignA += seq1[i - 1]
        alignB += '-'
        i -= 1
    
    while j > 0:
        alignA += '-'
        alignB += seq2[j - 1]
        j -= 1
    
    #Reversing the sequences as they were flipped from the traceback
    alignA = alignA[::-1]
    alignB = alignB[::-1]
    
    print("Finished traceback and sequences reversed.")
    matchString = ""
    for i in range (len(alignA)):
        if alignA[i] == alignB[i]:
            matchString += '|'
        elif alignA[i] != alignB[i]:
            if (alignA[i] == '-' or alignB == '-'):
                matchString += " "
            else:
                matchString += "*"
    alignScore = 0
    for i in range(len(matchString)):
        if matchString[i] == '|':
            alignScore += 1
        elif (matchString[i] == '*' or matchString[i] == " "):
            alignScore += -1
    
    return alignA, alignB, matchString, alignScore
output1, output2, matchString, alignScore, = NeedleWunsch(seq1, seq2, gap = -1, match = 1)

##Storing the outputs of seq1 and seq2 and other variables and printing the results
print("Alignment completed: ", output1[0:10] + '\n' + matchString[0:10]+ '\n' + output2[0:10])
print("Alignment Score: ", alignScore)

'''
##Adding the highlighting function. When there are two or more consercutive nucleotides that are not 
##matching between the two sequences, they will be represeneted in lowercase letters in the txt file.
##The txt file should be fluid and easy to read. It will have the seq1 and seq2 and its match string broken into chunks for easy reading.
##The difference will be in lowercase. And the sequence name will be stated in the beginnning of each sequence, even in the chunks
'''
def HighlightRegions(output1, output2):
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
print("Finished highlighting regions.")

def WriteToFile(output1, output2, matchString, seq1ID, seq2ID, output_file):
    print("Started writing alignment to text file.")
    #calls up highlightRegion function and stores the alignments in 
    highlightSeq1, highlightSeq2 = HighlightRegions(output1, output2)
    

    chunk_size = 60                 ##How wide the chunks of the alignment will be displayed as
    ##Had big issue with the match string going in between the sequence IDs, this is a good fix
    seqPadding = max(len(seq1ID), len(seq2ID)) + 2
    with open(output_file, "w") as file:             ##To make sure that the matchstring starts and is accurate to the alignment
        for start in range(0, len(output1), chunk_size):        ##Formatting how everything should be written into the file
            chunk1 = highlightSeq1[start:start + chunk_size]
            chunk2 = highlightSeq2[start:start + chunk_size]
            matchChunk = matchString[start:start + chunk_size]

            file.write(f"{seq1ID.ljust(seqPadding)}{chunk1}\n")
            file.write(f"{' ' * seqPadding}{matchChunk}\n")
            file.write(f"{seq2ID.ljust(seqPadding)}{chunk2}\n\n")

    print(f"Highlighted alignment written to {output_file}")
WriteToFile(output1, output2, matchString, seq1ID, seq2ID, output_file = "NW_Alignment.txt")

