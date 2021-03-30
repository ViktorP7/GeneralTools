#!/bin/bash

# Script to run kraken tool for ribosomal rDNA on all fastq files in a folder

# Requirements - Please install pigz, it is faster than gzip and gunzip on account of multithreading

# Command structure
# 16s_rDNA_KrakenRunner-01_11_20.sh FastqFileEnding PathToDB

# Define fastq file ending variable, and path to DB
FastqFileEnding=$1
PathToDB=$2


# Find the fastq files and get the amount of them
FastqFiles=(`ls | grep "$FastqFileEnding$"`)
NumFastqs=${#FastqFiles[@]}

# Initiate a for loop to loop thru the fastq files and process them 
for (( i=0; i<${NumFastqs}; i+=2)) 
do
	# Get each pair
	PairOne=${FastqFiles[$i]}
	PairTwo=${FastqFiles[$i+1]}

	# Create identifier prefix
	PairID=`echo $PairOne | awk '{ split($0, array, "_"); print array[1] }'`

	# Progress
	echo -e "\e Running sixess on $PairID... \e[0m"
	
	# Run sixess
	kraken2 --db $2 --paired $PairOne $PairTwo --report $PairID

	# Progress
	echo -e "\e[1;31m Pair $PairID fully processed, moving onto next... \e[0m"

done
	
	
