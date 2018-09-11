#!/bin/python3
# MEIN50010 2017 Assignment 2
# STUDENT ID 13460328
# Code to examine the genes (and some of their characteristics) in a given sequence

# Import sys module to get command line arguments
import sys

# Check right amount of arguments inputted
if len(sys.argv) < 6:
    print("Incorrect number of arguments, please try again using instructions below: \n"
          "Command line structure as follows: \n"
          "python mein50010-a2-13460328.py name.fasta name.gff table.txt output1.txt output2.txt value biopython \n"
          "name.fasta   The name of the fasta file you wish to examine \n"
          "name.gff     The annotation file in .gff format accompanying your fasta sequence \n"
          "table.txt    A table of the codon to amino acid translations \n"
          "output1.txt  The general output file, should be empty so you don't lose your own data \n"
          "output2.txt  A second output file for kmer related data \n"
          "value        An integer specifying the length of kmer you wish to examine. Optional, default is 3 \n"
          "biopython    Enter y/n or yes/no or 0/1 to confirm if you wish to use biopython or not. Optional, default no")
    quit()

# Define arguments as variables
seqFile = sys.argv[1]
annotationFile = sys.argv[2]
tableFile = sys.argv[3]
outputFile = sys.argv[4]
kcountFile = sys.argv[5]

# Check if input present for 6 and 7
if len(sys.argv) == 6:
    kinteger = int(3)
else:
    kinteger = int(sys.argv[6])

if len(sys.argv) == 6 or len(sys.argv) == 7:
    ynbiopython = "n"
else:
    ynbiopython = sys.argv[7]


# Find out if user wants to use biopython or my script to solve the problem
if ynbiopython == "n" or ynbiopython == "no" or ynbiopython == "0":

    # Define a function to process fasta file
    def get_seq(fastafile):

        # Open the sequence file for reading
        try:
            inputfastalines = open(fastafile, "r")
        # If that doesn't work, print below error message for file not found error
        except FileNotFoundError:
            print("could not open file {} for reading".format(fastafile))
            quit()

        # Initialise variable to store sequence
        seq = ""

        # Go thru each line of the fasta file and connect up the strings to form the whole sequence
        for line in inputfastalines:

            # Strip special characters from lines
            line = line.strip()

            # Skip the name, not interested in this case
            if ">" in line:
                continue

            # Add on the lines into the seq variable to get the whole seq
            else:
                seq = seq + str(line)

        # Close the sequence file, no longer needed
        inputfastalines.close()

        return seq

    # Call the above function to get sequence
    theSeq = get_seq(seqFile)

    # Function to calculate GC content of sequence
    def get_seq_gc(seq):

        # Calculate the proportion of GC in the sequence by counting up how much of each base is present
        a = seq.count("a") + seq.count("A")
        c = seq.count("c") + seq.count("C")
        g = seq.count("g") + seq.count("G")
        t = seq.count("t") + seq.count("T")

        # Need to add all up to get total count, add G and C for GC, and then divide to get the proportion
        totalbases = a + c + g + t
        gccontent = g + c
        proportiongc = float(gccontent) / float(totalbases)

        return proportiongc

    # Call above function to get GC content
    GCprop = get_seq_gc(theSeq)

    # Open file for writing
    fileOutput = open(outputFile, "w")

    # Write the GC proportion to file
    fileOutput.write("GC content: {:1.2f}\n\n".format(GCprop))

    # Define function to obtain a codon to amino acid translation table
    def get_translation_dict(codontable):

        # Open the table file for reading
        try:
            inputtable = open(codontable, "r")
        # If that doesn't work, print below error message for file not found error
        except FileNotFoundError:
            print("could not open file {} for reading".format(codontable))
            quit()

        # Make a dictionary to store codons and corresponding amino acids
        # Make lists to store both the amino acids and the bases
        tsldict = {}
        codons = []
        aminoacids = []

        # Go thru each line of the table and append the relevant things to the relevant lists
        for line in inputtable:

            line = line.strip()

            currentcodon = line.split("\t")[0]
            codons.append(currentcodon)

            currentacid = line.split("\t")[3]
            aminoacids.append(currentacid)

        # Close the codon table, no longer needed
        inputtable.close()

        # Pair both lists and enter them into the same dictionary
        for index in range(len(codons)):
            tsldict[codons[index]] = aminoacids[index]

        return tsldict

    # Call the dictionary function
    transDict = get_translation_dict(tableFile)

    # Lists to store gene names and gene seqs
    genenames = []
    geneseqs = []

    # Open the annotation file for reading
    try:
        inputAnnotation = open(annotationFile, "r")
    # If that doesn't work, print below error message for file not found error
    except FileNotFoundError:
        print("could not open file {} for reading".format(annotationFile))
        quit()

    # Must go thru the info for each gene
    for gene in inputAnnotation:

        gene = gene.strip()

        # To get gene name, chop the string at ID= and the semicolon ;, append gene names into a list
        geneName = gene.split("ID=")[1]
        geneName = geneName.split(";")[0]
        genenames.append(geneName)

        # Write current gene name to file
        fileOutput.write("Gene Name: {}\n".format(geneName))

        # Define function to get the current sequence
        def get_current_seq(gene, fullseq):

            # Get seq start and seq end positions for each gene from the gff file
            seqbits = gene.split("CDS")[1]
            seqbits = seqbits.split(".")[0]
            seqbits = seqbits.strip()
            seqstart = seqbits.split("\t")[0]
            seqend = seqbits.split("\t")[1]

            # Int the values, take away 1 from start value to account for python indexing, use to get relevant seq segment
            # Get the length of each gene
            # Append each gene into the seqeunces list
            seqstart = int(seqstart) - 1
            seqend = int(seqend)
            currseq = fullseq[seqstart:seqend]

            return(currseq)

        # Call function to get the current sequence
        currentSeq = get_current_seq(gene, theSeq)

        # Get length of the seuqnece, and append it into the sequence list
        currlen = len(currentSeq)
        geneseqs.append(currentSeq)

        # Write current sequence to file
        fileOutput.write("DNA sequence: {}\n".format(currentSeq))

        # Write current length to file
        fileOutput.write("Gene Length: {}\n".format(str(currlen)))

        # Get the GC content of each gene, similar as before
        currPropGC = get_seq_gc(currentSeq)

        # Write current GC content to file
        fileOutput.write("GC content: {:1.2f}\n".format(currPropGC))

        # Check which kind of strand the current gene is
        strand = gene.split("ID")[0]

        # Define a function to translate the current gene
        def get_strand_translation(currseq):

            currtsl = ""

            # Increment up in threes and thus check each codon against the dictionary and perform the translation
            increment = 3
            for base in range(0, len(currseq), increment):
                currcodon = (currseq[base:base + increment])
                currcodon = currcodon.upper()
                if currcodon in transDict:
                    currtsl = currtsl + transDict[currcodon]

            return (currtsl)

        if "+" in strand:

            # Call the translation function
            currentTrans = get_strand_translation(currentSeq)

            # Write translation to file
            fileOutput.write("Protein sequence: {}\n\n".format(currentTrans))

        elif "-" in strand:

            # Define a function to get the complimented sequence
            def get_compliment_seq(currseq):

                # First, reverse the currseq string, create a dictionary with the base compliments
                # Create a string variable to store the compliment
                # Go thru each base in the sequence and convert it to the compliment
                currseq = currseq[::-1]
                complimentbases = {"a":"t", "c":"g", "g":"c", "t":"a"}
                compliseq = ""
                for base in currseq:
                    if base in complimentbases:
                        compliseq = compliseq + complimentbases[base]

                return(compliseq)

            # Call function to get compliment sequence
            complimentSeq = get_compliment_seq(currentSeq)

            # Call function to get the translation done
            currentTrans = get_strand_translation(complimentSeq)

            # Write translation to file
            fileOutput.write("Protein sequence: {}\n\n".format(currentTrans))

    # Close the annotation file, not needed
    inputAnnotation.close()

    # Close the general output file
    fileOutput.close()

    # Dictionary to store total kmer counts
    ktotal = {}

    # Function to count kmers in a given sequence (from list of sequences)
    def countkmers(sequence, kinteger, kcounts):
        for index in range(len(sequence) - kinteger + 1):
            kmer = sequence[index:index + kinteger]
            if kmer in kcounts:
                kcounts[kmer] += 1
            else:
                kcounts[kmer] = 1
        return kcounts

    # Go thru each sequence in the geneseqs list and count up total kmers
    for sequence in range(len(geneseqs)):

        # Count the total
        countkmers(geneseqs[sequence], kinteger, ktotal)

    # Define a function to create a matrix filled with k counts
    def get_kmer_counts(genenames, ktotal, geneseqs, kinteger):

        # Create a matrix that corresponds to the amount of gene names and total kmers present
        kmermatrix = [["" for x in range(len(genenames) + 1)] for y in range(len(ktotal) + 1)]

        # Fill the matrix with the gene names and also "kmer" in top left
        kmermatrix[0][0] = "kmer"

        for index in range(len(genenames)):
            kmermatrix[0][index + 1] = genenames[index]

        # Get names of ktotal, place in a list, add to matrix
        knames = ktotal.keys()

        for index in range(len(knames)):
            kmermatrix[index + 1][0] = knames[index]

        # Go thru each sequence in the geneseqs list and count up individual kmers and insert into the matrix
        for sequence in range(len(geneseqs)):

            # Create a dictionary to store kmer counts
            kcounts = {}

            # Count for each individual sequence
            countkmers(geneseqs[sequence], kinteger, kcounts)

            # Go thru each kmer in the corresponding matrix slot, if present in kcounts, transfer the value over
            # If not in kcounts, note down a zero for that spot
            for index in range(len(knames)):
                if kmermatrix[index + 1][0] in kcounts:
                    kmermatrix[index + 1][sequence + 1] = kcounts[kmermatrix[index + 1][0]]
                else:
                    kmermatrix[index + 1][sequence + 1] = 0

        return kmermatrix

    # Call the function and create a k count matrix
    kMatrix = get_kmer_counts(genenames, ktotal, geneseqs, kinteger)

    # Open k count file for writing
    filekCount = open(kcountFile, "w")

    # Write the k count matrix to file
    for line in kMatrix:
        filekCount.write("\n")
        for thing in line:
            filekCount.write("%8s" %(thing))

    # Close the k count file
    filekCount.close()

    # Notify end of program
    print("All tasks successfully completed, please check your output files using Notepad++")

# If biopython is yes, then use biopython module to solve issue
elif ynbiopython == "y" or ynbiopython == "yes" or ynbiopython == "1":
    print("You have elected to use the biopython module")

else:
    print("Input not recognised, please carefully read input parameters below and try again: \n"
          "Command line structure as follows: \n"
          "python mein50010-a2-13460328.py name.fasta name.gff table.txt output1.txt output2.txt value biopython \n"
          "name.fasta   The name of the fasta file you wish to examine \n"
          "name.gff     The annotation file in .gff format accompanying your fasta sequence \n"
          "table.txt    A table of the codon to amino acid translations \n"
          "output1.txt  The general output file, should be empty so you don't lose your own data \n"
          "output2.txt  A second output file for kmer related data \n"
          "value        An integer specifying the length of kmer you wish to examine \n"
          "biopython    Enter y/n or yes/no to confirm if you wish to use biopython or not")
    quit()
