# Tool to parse vcf files and produce a merged vcf table
# Author Viktor Perets
# Version 1.0 15-09-21
# Must specify filters as based on initial program parameters

###################
#### FUNCTIONS ####
###################

# Define function to read in contents of all files
def get_contents(fileEnding, nameOfDir):

	# Create dictionary and arrays to store file contents and names
	contentArray = []
	nameArray = []
	fileDict = {}

	# Loop thru items in dir
	for item in os.listdir(nameOfDir): 

		# Check for "csq.vcf" extension
		if item.endswith(fileEnding): 

			# Open file for reading
			file = open(item, "r")

			# Read in the content
			content = file.readlines()

			# Close current file
			file.close()
			
			# Store contents
			contentArray.append(content)
			
			# Store name
			nameArray.append(item)
			
	# Pair both lists and enter them into the same dictionary
	for index in range(len(nameArray)):
		fileDict[nameArray[index]] = contentArray[index]
			
	return(fileDict)

# Define function to store positions of each line in the file
def get_all_positions(files):
	
	# Arrays and dictionaries to store all positions, lines and file names
	filePosLine = {}
	positionsarray = []
	linearray = []
	
	# Loop thru each file in all file
	for file in files:
		
		# Access via key
		currFile = files[file]
	
		# Loop thru each line of the file
		for line in currFile:
		
			# Strip the line
			line = line.strip()
		
			# Skip if hash
			if line[0] == "#":
			
				continue
		
			# Otherwise store the position
			else:
			
				# Store split line
				splitLine = line.split("\t")
				
				# Store the allele position
				position = "{}:{}".format(file, splitLine[1])
				
				# Store position in the array
				positionsarray.append(position)
				
				# Store the line
				linearray.append(line)
				
		# Pair both lists and enter them into the same dictionary
		for index in range(len(positionsarray)):
			filePosLine[positionsarray[index]] = linearray[index]
						
					
	return(filePosLine)

# Define function to store positions of each line in the file
def get_unique_positions(files):
	
	# Array to store unique positions
	positionsarray = []

	# Loop thru each file in all file
	for file in files:
		
		# Access via key
		currFile = files[file]
	
		# Loop thru each line of the file
		for line in currFile:
		
			# Strip the line
			line = line.strip()
		
			# Skip if hash
			if line[0] == "#":
			
				continue
		
			# Otherwise store the position
			else:
			
				# Store split line
				splitLine = line.split("\t")
				
				# Store the allele position
				position = splitLine[1]
				
				# Check if present in array
				if position in positionsarray:
					
					continue
				
				else:
				
					# Store position in the array
					positionsarray.append(position)
						
					
	return(positionsarray)

# Define a function to create a matrix filled with the positions and note presence in each file
def get_position_presence_matrix(files, positions, uniquePositions):

	# Create a matrix that corresponds to the amount of positions and files present
	csqMatrix = [["NA" for x in range(len(files) + 1)] for y in range(len(uniquePositions) + 1)]

	# Fill the matrix with the file names and also "Position" in top left
	csqMatrix[0][0] = "Position"
	
	# Get names of the files
	fnames = list(files.keys())

	for index in range(len(fnames)):
		csqMatrix[0][index + 1] = fnames[index]

	# Fill matrix with the positions
	for index in range(len(uniquePositions)):
		csqMatrix[index + 1][0] = uniquePositions[index]

	# Go thru each position and check if present in each file
	for position in positions:
	
		# Initialise variables to store coordinates of a position
		coordsX = None
		coordsY = None
					
		# Store split key
		splitPos = position.split(":")
		fileName = splitPos[0]
		posName = splitPos[1]
		
		# Find the index of the file and position in the respective arrays
		coordsY = 1 + fnames.index(fileName)
		coordsX = 1 + uniquePositions.index(posName)
		
		# Now access the line stored under the key and split by tabs
		splitLine = positions[position].split("\t")
		
		# Store the quality and csq info
		info = splitLine[7]
					
		# Initialise DP4, DP, CSQ and MQ variables
		DP4 = None
					
		DP = None
					
		CSQ = None
		
		MQ = None
					
		# Initialise a status (set to pass initially)
		status = "PASS"
					
		# Store the DP, DP4, BCSQ and MQ quality parameters
		for score in info.split(";") :
						
			if "DP4" in score :
							
				DP4 = score
				
			elif "BCSQ" in score:
						
				CSQ = score
						
			elif "DP" in score :
							
				DP = score
			
			elif "MQB" in score:
			
				continue
				
			elif "MQSB" in score:
			
				continue
			
			elif "MQ0F" in score:
			
				continue
			
			elif "MQ" in score:
			
				MQ = score
							
			
							
		# Get DP4 values by splitting 
		DP4vals = DP4.split("=")[1]
					
		DP4Indivals = DP4vals.split(",")
					
		# Initialise variables for ref and alt forward and reverse, and variables for the total and proportions
		refFw = float(DP4Indivals[0])
		refRev = float(DP4Indivals[1])
		altFw = float(DP4Indivals[2])
		altRev = float(DP4Indivals[3])
		DP4total = refFw + refRev + altFw + altRev
		refProp = float((refFw + refRev) / DP4total) 
		altProp = float((altFw + altRev) / DP4total) 
					
		# Fail the SNP if DP is below 30
		if float(DP.split("=")[1]) < 10 :
					
			status = "FAIL-DP"
			
		# Fail the MQ if it's below 30
		if float(MQ.split("=")[1]) < 20:
		
			status = "FAIL-MQ"
					
		# If DP4 alt forward and reverse are each <4, fail
		elif altFw < 2 or altRev <2:
						
			status = "FAIL-Support"
						
		# If the proportion of DP4 < 0.95, fail
		elif altProp < 0.95:
					
			status = "FAIL-SupportProp"
					
		# Store it as a string and place inside matrix
		quality = "{}:{}".format(CSQ, status)
					
		csqMatrix[coordsX][coordsY] = quality

	return csqMatrix
###################
## END FUNCTIONS ##
###################
	
# Import modules
import os, sys, datetime

# Check right amount of arguments inputted
if len(sys.argv) < 1:
	print("Incorrect number of arguments, please try again using instructions below: \n"
		  "Command line structure as follows: \n"
		  "python Get-VCF-Positions.py directory_path \n"
		  "directory_path   The full path of the directory your files are in, or type cwd for current"
		  )
	quit()

# Define arguments as variables
dirpath = sys.argv[1]

# Check what is in dirpath and set directory accordingly
if dirpath == "cwd":

	dir_name = os.getcwd()

else:
	
	dir_name = dirpath

# Set Directory and Extension names
extension = ".vcf"

# Change directory from working dir to dir with files
os.chdir(dir_name) 

# Get the file contents
allFiles = get_contents(extension, dir_name)

# Get the positions
allPositions = get_all_positions(allFiles)

# Get the unique positions
uniques = get_unique_positions(allFiles)

# Create empty arrays to store sorted unique positions
intUniques = [None] * len(uniques)
sortedUniques = [None] * len(uniques)

# Make uniques into ints and store
for index in range(len(uniques)):
	intUniques[index] = int(uniques[index])
	
# Sort the ints
intUniques.sort()

# Convert back to strings and place in sorted uniques
for index in range(len(intUniques)):
	sortedUniques[index] = str(intUniques[index])
	
# Get the matrix
bigMatrix = get_position_presence_matrix(allFiles, allPositions, sortedUniques)

# Create an output file
outputFile = "ComVCF_{}.tsv".format(datetime.datetime.today().strftime('%Y-%m-%d'))
print(outputFile)

# Open output for writing
output = open(outputFile, "w")

# Write the matrix to file
for line in bigMatrix:
	
	# New line for each line
	output.write("\n")
	
	# For each thing in each line
	for thing in line:
	
		# Separate by commas to make .csv file
		output.write("{}\t".format(thing))
	
# Close the output file
output.close()
				
# Print completion
print("Task complete, please check output file in your current directory for results")					

			
				
			
			
		
		