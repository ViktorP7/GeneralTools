# Import modules
import sys

# Check right amount of arguments inputted
if len(sys.argv) < 2:
    print("Incorrect number of arguments, please try again using instructions below: \n"
          "Command line structure as follows: \n"
          "python CompareConsequences20-08-18.py File1 File2 \n"
          "File1 The first consequence file you want to compare \n"
		  "File2 The second consequence file you want to compare"
          )
    quit()

# Define arguments as variables
file1 = sys.argv[1]
file2 = sys.argv[2]

# Create an output file
outputFile = "CCsq.txt"

# Open output file for writing
output = open(outputFile, "w")

def find_difference(File1, File2):

	# Open file 1 for reading
	with open(File1) as f1:
	
		# Open file 2 for reading
		with open(File2) as f2:
		
			# Write out the input file name 
			output.write("#Source:{}\n".format(File1))
		
			# Loop thru the lines of file 1
			for line1 in f1:
		
				# Skip line if there is a hash at the start of the line
				if line1[0] == "#" :
			
					continue
			
				# Split the line by tabs
				else:
				
					# Store split line
					splitLine1 = line1.split("\t")
				
					# Store the allele position
					position1 = splitLine1[1]
				
					# Loop thru f2
					for line2 in f2:
							
							# Skip line if there is a hash at the start of the line
							if line2[0] == "#" :
			
								continue
			
							# Split the line by tabs
							else:
				
								# Store split line
								splitLine2 = line2.split("\t")
				
								# Store the allele position
								position2 = splitLine2[1]
								
								# Check if the positions equal
								if position1 != position2:
				
									# Write full line1 to output 
									output.write("{}\n".format(line1))
									
									
# Find what is in 1 and not in 2
find_difference(file1, file2)

# Find what is in 2 and not in 1
find_difference(file2, file1)
		
			
					
# Print completion
print("Task complete, please check output file in your current directory for results")					
# Close the output file
output.close()
			
				
			
			
		
		