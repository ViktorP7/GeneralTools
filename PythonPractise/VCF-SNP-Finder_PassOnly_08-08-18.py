# Import modules
import os, gzip, sys

# Check right amount of arguments inputted
if len(sys.argv) < 1:
    print("Incorrect number of arguments, please try again using instructions below: \n"
          "Command line structure as follows: \n"
          "python directory_path \n"
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
extension = ".gz"

# Create an output file
outputFile = "SNPFinder.txt"

os.chdir(dir_name) # change directory from working dir to dir with files

# Open output file for writing
output = open(outputFile, "w")

for item in os.listdir(dir_name): # loop thru items
    
	if item.endswith(extension): # check for ".gz" extension
	
		# Initiate a SNP counter
		SNPcounter = 0
		
		# Write out the input file name 
		output.write(">Source:{}\n#Count\tPosition\tReference\tAlternate\tDP4\tDP\tStatus\n".format(item))
		
		# Unzip file
		file = gzip.open(item, "rb")
		
		# Read in the content
		content = file.readlines()
		
		# Close current file
		file.close()
		
		# Loop thru the content 
		for line in content :
		
			line = line.strip()
		
			# Skip line if there is a hash at the start of the line
			if line[0] == "#" :
			
				continue
			
			# Split the line by tabs
			else:
				
				# Store split line
				splitLine = line.split("\t")
				
				# Store the allele position
				position = splitLine[1]
				
				# Store the reference allele
				refallele = splitLine[3]
				
				# Store the alternative allele
				altallele = splitLine[4]
				
				# Store the quality information
				quality = splitLine[7]
				
				# Check if the ref and alt are identical
				if altallele == "." :
					
					continue
					
				else :
				
					# Initialise a status (set to pass initially)
					status = "PASS"
					
					# Initialise DP4 and DP variable
					DP4 = None
					
					DP = None
					
					# Store the DP and DP4 quality parameters
					for score in quality.split(";") :
						
						if "DP4" in score :
							
							DP4 = score
						
						elif "DP" in score :
							
							DP = score
					
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
					if float(DP.split("=")[1]) < 30 :
					
						status = "FAIL"
					
					# If DP4 alt forward and reverse are each <4, fail
					elif altFw < 4 or altRev <4:
						
						status = "FAIL"
						
					# If the proportion of DP4 < 0.95, fail
					elif altProp < 0.95:
					
						status = "FAIL*"
						
					# Write only if pass
					if status == "PASS":
					
						# Update SNP counter
						SNPcounter = SNPcounter + 1
					
						# Write out the number of current SNP position, reference alt allele, quality and status if the ref and alt differ
						output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(SNPcounter, position, refallele, altallele, DP4, DP, status))


# Print completion
print("Task complete, please check output file in your current directory for results")					
# Close the output file
output.close()
			
				
			
			
		
		