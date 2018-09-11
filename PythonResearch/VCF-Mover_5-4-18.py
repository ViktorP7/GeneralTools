#!/bin/python3

# Import sys module to get command line arguments
import sys

# Check right amount of arguments inputted
if len(sys.argv) < 1:
    print("Incorrect number of arguments, please try again using instructions below: \n"
          "Command line structure as follows: \n"
          "python vcfmover5_4_18.py dump.names \n"
          "dump.names Names of VCF files you wish to be moved and dumped into a dump folder")
    quit()
	
# Define arguments as variables
dumpFile = sys.argv[1]

# Define a function to process dump file
def get_dump(dumpfile):

    # Open the dump file for reading
    try:
        inputdumplines = open(dumpfile, "r")
    # If that doesn't work, print below error message for file not found error
    except FileNotFoundError:
        print("could not open file {} for reading".format(dumpfile))
        quit()

    # Initialise array to store dump
    dump = []

    # Go thru each line of the dump file
    for line in inputdumplines:

        # Strip special characters from lines
        line = line.strip()
        # Add on the lines into the dumpo
        dump.append(line)

    # Close the dump file, no longer needed
    inputdumplines.close()

    return dump

# Call the above function to get sequence
theDump = get_dump(dumpFile)

# Take away the header from the array
del theDump[0]

# New array to store dump names
realDump = []

# For each element in the dump array
for element in theDump:

	newelement = element.split(" ")[1]
	
	realDump.append(newelement)

# Import os 
import os

# Create a folder to store the removed files in current directory
os.makedirs("Removed_files")

# Import pathlib and shutil
from pathlib import Path
import shutil

# Loop thru each file name in the real dump array
for name in realDump:

	if os.path.isfile(name):
		
		shutil.move(name, "Removed_files" + "/" + name)
		
		

	