# Script to extract gene info for SNPs outputted by the SNP-Finder tool
# Tailored to with MAP K10 annotations only for now

# Import modules
import os, sys

# Check right amount of arguments inputted
if len(sys.argv) < 3:
	print("Incorrect number of arguments, please try again using instructions below: \n"
		  "Command line structure as follows: \n"
		  "python directory_path SNPFinderFile AnnotationFile \n"
		  "directory_path   The full path of the directory your files are in, or type cwd for current \n"
		  "SNPFinderFile	The output file from the VCF-SNP-Finder tool \n"
		  "AnnotationFile	The annotation file in .gff or .gff3 format for the reference genome"
		  )
	quit()

# Define arguments as variables
dirpath = sys.argv[1]
SNPFile = sys.argv[2]
AnnoFile = sys.argv[3]

# Check what is in dirpath and set directory accordingly
if dirpath == "cwd":

	dir_name = os.getcwd()

else:
	
	dir_name = dirpath

# Create an output file
outputFile = "GeneFinder.txt"

os.chdir(dir_name) # change directory from working dir to dir with files

# Read in the SNP lines
SNPLines = open(SNPFile, "r")

# Read in the .gff or .gff3 files
AnnoImp = open(AnnoFile, "r")
AnnoLines = AnnoImp.readlines()

# Open output file for writing
#output = open(outputFile, "w")

	
# Initialise variables to store current name and position
name = ""
position = ""

# Go thru each line of SNP file
for line in SNPLines:

	# Strip special characters from lines
	line = line.strip()

	# Name the name
	if ">" in line:

		name = line
		print(name)
		
		# Write the name to output
		#output.write("{}\n".format(name))
		
	# Skip the info column
	elif "#" in line:
	
		continue
		

	# Get the position by chopping
	else:

		splitter = line.split("\t")
		
		# Float the position
		position = float(splitter[1])
			
		# Loop thru the annotations file
		for anno in AnnoLines:
			
			# Skip the header line
			if anno[0] == "#":
					
				continue
			
			# Skip line corresponding to entire genome
			elif "region" in anno:
				
				continue
			
			# Skip anything that doesn't include RefSeq - RefSeq and Protein Homology lines will not be skipped
			elif "RefSeq" not in anno:
				
				continue
					
			else:
				
				# Split the annotation line by tabs
				choppedAnno = anno.split("\t")
				
				# Check if the position is within the start and end position for the current gene
				if position >= float(choppedAnno[3]) and position <= float(choppedAnno[4]):
					
					# If the line is gene
					if choppedAnno[2] == "gene":
						
						# Split the info column
						info = choppedAnno[8].split(";")
						
						# Store the name of the gene
						print(info[2])
					
					# If it is a protein line
					elif choppedAnno[1] == "Protein Homology":
						
						# Split info column
						info = choppedAnno[8].split(";")
						
						# Store the info about the protein product
						print(info[6])
		


# Close the files
SNPLines.close()
AnnoImp.close()
