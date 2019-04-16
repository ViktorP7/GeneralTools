#!/bin/python3

# Import sys module to get command line arguments sorted
import sys

# Check right amount of arguments inputted (need two)
if len(sys.argv) != 2:
    print("incorrect number of arguments")
    quit()

# Define arguments as variables
fileName = sys.argv[1]

# Open the input file for reading
try:
    inputFastaLines = open(fileName, "r")
# If that doesn't work, print below error message for file not found error
except FileNotFoundError:
    print("could not open file {} for reading".format(fileName))
    quit()

# Create output file	
fastaOutput = "Chopped" + fileName

# Open log file for writing
outputFasta = open(fastaOutput, "w")

# Initialise arrays to store sequence name and array for sequences
seqname = []
seqs = []

# Store, count up and write out the amount of sequences present in the input file
someseq = ""
for line in inputFastaLines:

    # Remove new line character from line
    line = line.strip()

    if ">" in line:

        # Check whether there is a sequence that needs storing
        if someseq != "":
            seqs.append(someseq)

        # Append current line as a string to the seqname array, reset someseq variable
        seqname.append(str(line))
        someseq = ""

    # Concatenate the current line (string) to someseq
    else:
        someseq = someseq + str(line)

# Store the last sequence
seqs.append(someseq)

# Close input file as it is no longer needed
inputFastaLines.close()

# Reuse the names and sequences arrays from before
# Chop off the excess info from names (so anything after a "_")
# Initialise new name array
seqnamenew = []
for names in seqname:
    seqnamenew.append(names.split("_")[0])

# Go thru each sequence and write it to file
for line in range(len(seqnamenew)):
    outputFasta.write("{}\n".format(seqnamenew[line]))
    currentseq = seqs[line]
    outputFasta.write(currentseq + "\n")


# Close the reformatted fasta file
outputFasta.close()

# Print completion message
print("Status: All tasks completed! Please check working directory for desired chopped fasta!")