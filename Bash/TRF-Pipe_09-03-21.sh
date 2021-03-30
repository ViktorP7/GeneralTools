#!/bin/bash

# Pipeline for using Tandem Repeats Finder (TRF) on assmbled WGS data
# Pipe MUST be run from DESKTOP, else spades will not work properly!!! 
# Designed around running from Desktop
# Command structure
# TRF-Pipe_09-03-21.sh FastqFileEnding PathToFastqs

# Define fastq file ending variable, path to reference and annotations
FastqFileEnding=$1
PathToFastqs=$2

# Find the fastq files and get the amount of them
FastqFiles=(`ls $PathToFastqs | grep "$FastqFileEnding$"`)
NumFastqs=${#FastqFiles[@]}

# Make a directory to house the TRF results
mkdir UbuntuSharedFolder/ResultsTRF

# Initiate a for loop to loop thru the fastq files and process them 
for (( i=0; i<${NumFastqs}; i+=2)) 
do
	# Get each pair
	PairOne=${FastqFiles[$i]}
	PairTwo=${FastqFiles[$i+1]}

	# Create identifier prefix
	PairID=`echo $PairOne | awk '{ split($0, array, "_"); print array[1] }'`

	# Progress
	echo -e "Copying files to Desktop and beginning Spades assembly... \e[0m"
	
	# Copy files to Desktop
	cp $PathToFastqs$PairID* . 
	
	# Make assembly path name
	AssemblyPath=$PairID"_assembly"
	
	# Run Spades
	spades.py --pe1-1 $PairOne --pe1-2 $PairTwo -o $AssemblyPath

	# Variable to store scaffolds
	FastaFile=$PairID".fa"
	
	# Concatenate the contigs
	cat $AssemblyPath"/scaffolds.fasta" | sed -e '1!{/^>.*/d;}' | sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{80\}/&\n/g' > $FastaFile

	# Progress
	echo -e "\e[1;31m Assembly complete, running TRF... \e[0m"
	
	# Run TRF
	trf $FastaFile 2 5 7 80 10 50 2000 -l 10 -h
	
	# Create an output for the TRF
	FinalOut=$PairID".dat"
	
	# Format the TRF output
	tail -n +16 $FastaFile".2.5.7.80.10.50.2000.dat" > $FinalOut
	
	# Move final output to new folder
	mv $FinalOut UbuntuSharedFolder/ResultsTRF

	# Progress
	echo -e "\e[1;31m Removing intermediate files... \e[0m"

	# Remove unecessary files
	rm -rf $AssemblyPath
	rm $PairID*

	# Progress
	echo -e "\e[1;31m Pair $PairID fully processed, moving onto next... \e[0m"

done
