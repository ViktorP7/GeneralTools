#!/bin/python3
# MEIN50010 2017 Assignment 1
# STUDENT ID 13460328
# Code to analyse and reformat a FASTA sequence
# I've left in print statements throughout the code for witnessing progress, not sure if you wanted them removed or not

# Import sys module to get command line arguments sorted
import sys

# Check right amount of arguments inputted (need two)
if len(sys.argv) != 3:
    print("incorrect number of arguments")
    quit()

# Define arguments as variables
fileName = sys.argv[1]
positiveInteger = float(sys.argv[2])

# Create names for output files by chopping off ".txt" or ".fa" from file name and replacing with appropriate
if ".txt" in fileName:
    statsOutput = fileName.replace(".txt", ".log")
    fastaOutput = fileName.replace(".txt", ".out")
elif ".fa" in fileName:
    statsOutput = fileName.replace(".fa", ".log")
    fastaOutput = fileName.replace(".fa", ".out")
else:
    print("Unknown file format, please submit as .txt or .fa")
    quit()

# Open the input file for reading
try:
    inputFastaLines = open(fileName, "r")
# If that doesn't work, print below error message for file not found error
except FileNotFoundError:
    print("could not open file {} for reading".format(fileName))
    quit()

# Open log file for writing
outputStats = open(statsOutput, "w")

# Write out the input file name at the start of the new log file
outputStats.write("Source: {}\n".format(fileName))

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

# Print information about sequences stored to screen
print("Number of sequence names: " + str(len(seqname)))
print("Number of sequences: " + str(len(seqs)))

# Write number of sequences to log file
outputStats.write("Number of sequences: {}\n".format(len(seqs)))

# Close input file as it is no longer needed
inputFastaLines.close()


# Define min and max functions
def min(a, b):
    if a < b:
        return a
    else:
        return b


def max(a, b):
    if a > b:
        return a
    else:
        return b


# Define min and max variables, create array to store lengths of sequences
myMax = 0
myMin = 9999999999999999999999999999999
seqlengths = []

# Use above functions to calculate min and max values for the sequences, ignoring spaces
# Also append sequence lengths into the seqlengths array
for sequence in seqs:
    myLen = len(sequence) - sequence.count("-")
    myMin = min(myMin, myLen)
    myMax = max(myMax, myLen)
    seqlengths.append(myLen)

# Divide the sum of the sequence length to get the average value
averageseq = sum(seqlengths) / len(seqlengths)

# Print to screen
print("Min is: " + str(myMin))
print("Max is: " + str(myMax))

# Print the average value and sequence lengths vector to screen
print("Average is: " + str(averageseq))

# Write min and max and average length to log file
outputStats.write("Shortest: {}\n".format(myMin))
outputStats.write("Longest: {}\n".format(myMax))
outputStats.write("Average: {:1.1f}\n".format(averageseq))


# Assign the length of the first sequence to a variable as the default alignment length
alnLen = len(seqs[0])

# If the alignment is not the same for all sequences, -1 is assigned to alnLen to show no alignment and loop breaks
# Else, assuming everything is the same, the alignment length remains as is
for sequence in seqs:
    if alnLen != len(sequence):
        alnLen = -1
        break

# Print alignment length to screen
print("Alignment length is: " + str(alnLen))

# Store alignment length in the log output file
outputStats.write("Aligned: {}\n".format(alnLen))

# Reuse the names and sequences arrays from before
# Chop off the excess info from names (so anything after a blank space)
# Initialise new name array
seqnamenew = []
for names in seqname:
    seqnamenew.append(names.split(" ")[0])

# Get rid of dashes in sequence
# Initialise new array for new sequences
newseqs = []
for sequence in seqs:
    seqsplit = (sequence.split("-"))
    newseqs.append("".join(seqsplit))

# Initialise an empty variable type to store type of sequence
type = ""

# Check using first sequence(should apply to all anyway), if it is DNA, RNA or PROTEIN
# Note that all of the bases below must be present in the sequence for it to be classified 
# The RNA or DNA is checked not to contain a random amino acid, in this case "W",
# to make sure a protein containing the letters for DNA/RNA bases is not mislabelled
if ("A" in newseqs[0]) and ("T" in newseqs[0]) and ("G" in newseqs[0]) and ("C" in newseqs[0]) and (
    "W" not in newseqs[0]):
    type = "DNA"
    print("The type is: " + type)

# Check if it might be RNA
elif ("A" in newseqs[0]) and ("U" in newseqs[0]) and ("G" in newseqs[0]) and ("C" in newseqs[0]) and (
    "W" not in newseqs[0]):
    type = "RNA"
    print("The type is: " + type)

# If it's not RNA or DNA, it must be a protein
else:
    type = "PROTEIN"
    print("The type is: " + type)

# Write to log file the type
outputStats.write("Sequence Type: {}\n".format(type))

# Create dictionaries to correspond to each base or amino acid
basesRNA = {}
basesDNA = {}
aminoAcids = {}

# Refer back to type from earlier to determine which count is more appropriate
# Going thru each sequence in the sequences array, and thru each base
# Of every sequence, store the base in the dictionary upon first time
# Appearing, and afterwards add 1 to its count
# Sort the dictionary in alphabetical order of keys
# Create array to store the counts of each key
# Verify that correct bases were counted for RNA and DNA
# Append the key values to array, get the sum total
# Calculate proportion by dividing each key by the total of residues/bases present
# Print info to screen and store in log file
if type == "RNA":
    for sequence in newseqs:
        for base in range(len(sequence)):
            if sequence[base] in basesRNA:
                basesRNA[sequence[base]] += 1
            else:
                basesRNA[sequence[base]] = 1
    sortedBases = sorted(basesRNA.keys())
    rnaCounts = []
    for key in sortedBases:
        rnaCounts.append(basesRNA.get(key))

    rnatotal = sum(rnaCounts)
    rnaprop = []
    for value in rnaCounts:
        rnaprop.append(float(value) / float(rnatotal))
    print("RNA count is: " + str(rnaCounts))
    print("RNA total is: " + str(rnatotal))
    print("RNA proportion is: " + str(rnaprop))
    outputStats.write("Distribution: ")
    for prop in rnaprop:
        outputStats.write("{:1.2f} ".format(prop))
    outputStats.write("\n")




elif type == "DNA":
    for sequence in newseqs:
        for base in range(len(sequence)):
            if sequence[base] in basesDNA:
                basesDNA[sequence[base]] += 1
            else:
                basesDNA[sequence[base]] = 1

    try:
        basesDNA.pop("Y") # Get rid of annoying thing that shouldn't be here
    except KeyError:
        pass


    sortedBases = sorted(basesDNA.keys())
    dnaCounts = []
    for key in sortedBases:
        dnaCounts.append(basesDNA.get(key))

    dnatotal = sum(dnaCounts)
    dnaprop = []
    for value in dnaCounts:
        dnaprop.append(float(value) / float(dnatotal))
    print("DNA count is: " + str(dnaCounts))
    print("DNA total is: " + str(dnatotal))
    print("DNA proportion is: " + str(dnaprop))
    outputStats.write("Distribution: ")
    for prop in dnaprop:
        outputStats.write("{:1.2f} ".format(prop))
    outputStats.write("\n")


elif type == "PROTEIN":
    for sequence in newseqs:
        for base in range(len(sequence)):
            if sequence[base] in aminoAcids:
                aminoAcids[sequence[base]] += 1
            else:
                aminoAcids[sequence[base]] = 1

    sortedAcids = sorted(aminoAcids.keys())
    acidCounts = []
    for key in sortedAcids:
        acidCounts.append(aminoAcids.get(key))

    acidtotal = sum(acidCounts)
    acidprop = []
    for value in acidCounts:
        acidprop.append(float(value) / float(acidtotal))
    print("Acid count is: " + str(acidCounts))
    print("Acid total is: " + str(acidtotal))
    print("Acid proportion is: " + str(acidprop))
    outputStats.write("Distribution: ")
    for prop in acidprop:
        outputStats.write("{:1.2f} ".format(prop))
    outputStats.write("\n")

# Try find largest common substring across all sequences
# Note that this code will count spaces ("-") as a common thing
# First, check if the sequences are actually aligned, if not aligned, this operation cannot be performed
# If aligned, initialise an empty variable to store the common sequence
# Compare all sequences in the array and try find the common sequence, store in the empty variable
if alnLen == -1:
    print("Sequences out of alignment, can't find common substring")
    outputStats.write("Common: Sequences out of alignment, can't find common substring\n")
else:
    def commonstring(seqs):
        common = ""
        if len(seqs) > 1 and len(seqs[0]) > 0:
            for sequence in range(len(seqs[0])):
                for comparison in range(len(seqs[0]) - sequence + 1):
                    if comparison > len(common) and is_a_substring(seqs[0][sequence:sequence + comparison], seqs):
                        common = seqs[0][sequence:sequence + comparison]
        return common


    def is_a_substring(find, seqs):
        if len(seqs) < 1 and len(find) < 1:
            return False
        for sequence in range(len(seqs)):
            if find not in seqs[sequence]:
                return False
        return True


    print("Common substring is: " + str(commonstring(seqs)))
    outputStats.write("Common: {}\n".format(str(commonstring(seqs))))

# Close the log file, it is done now
outputStats.close()

# To write the reformatted fasta sequences to a new file, open the file for writing
outputFasta = open(fastaOutput, "w")

# Store desired line length value here
linevalue = int(positiveInteger)

# Go thru each sequence and write it to file
for line in range(len(seqnamenew)):
    outputFasta.write("{}\n".format(seqnamenew[line]))
    currentseq = newseqs[line]
    while len(currentseq) > 0:
        outputFasta.write(currentseq[:linevalue] + "\n")
        currentseq = currentseq[linevalue:]


# Close the reformatted fasta file
outputFasta.close()

# Print completion message
print("Status: All tasks completed! Please check working directory for desired .log and .out files!")