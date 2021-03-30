#!/bin/bash

# Alignment and consequence calling of paired raw reads
# Contains adapted code from ProcessRawReads_09-08-18 by Joe Crispell
# Author - Viktor Perets 16/08/18
# For bcftools csq to work properly, your fasta and annotations should be located in the working directory
# State exact names of reference and gff, not "."
# 24-04-19 Updated to remove unecessary intermediates, speeded up unzipping
# 29-10-20 Removed trimming for the moment, updated parameters
# 25-11-20 Changed aligner to use unpaired data (as well as other steps) - CIT RNA seqs are unpaired replicates!!!

# Requirements - Please install pigz, it is faster than gzip and gunzip on account of multithreading

# Command structure
# ConsequencePipe_16-08-18.sh FastqFileEnding PathToReference PathToAnnotations

# Define fastq file ending variable, path to reference and annotations
FastqFileEnding=$1
PathToReference=$2
PathToAnnotations=$3

# Find the fastq files and get the amount of them
FastqFiles=(`ls | grep "$FastqFileEnding$"`)
NumFastqs=${#FastqFiles[@]}

# Initiate a for loop to loop thru the fastq files and process them 
for (( i=0; i<${NumFastqs}; i+=1)) 
do
	# Get each pair
	PairOne=${FastqFiles[$i]}

	# Create identifier prefix
	PairID=`echo $PairOne | awk '{ split($0, array, "_"); print array[1] }'`

	# Progress
	echo -e "Beginning alignment of trimmed files... \e[0m"

	# Variables to store sam, bam sorted bam and indexed bam
	Samfile=$PairID"_aln-se.sam"
	Bamfile=$PairID".bam"
	SortedBam=$PairID"_srtd.bam"
	IndexedBam=$SortedBam".bai"
	NDupBam=$PairID"_srtd_ndup.bam"

	# Align using bwa aln
	FILE1SAI=$PairOne".sai"
	bwa aln -t 4 $PathToReference $PairOne > $FILE1SAI

	# Generate SAM file
	# -f Flag: SAM file to output results
	bwa samse -f $Samfile $PathToReference $FILE1SAI $PairOne 
	
	# Remove unwanted files
	rm $FILE1SAI
	
	samtools view -b $Samfile > $Bamfile
	samtools sort $Bamfile -o $SortedBam
	samtools index $SortedBam
	samtools rmdup $SortedBam $NDupBam

	# Progress
	echo -e "\e[1;31m Alignment complete, creating bcf file... \e[0m"

	# Variable to store bcf file
	BcfFile=$PairID".bcf"

	# Create a bcf file with samtools
	samtools mpileup --min-MQ 10 --uncompressed --fasta-ref $PathToReference $NDupBam > $BcfFile

	# Progress
	echo -e "\e[1;31m Bcf file created, calling variants... \e[0m"

	# Initialise variable to store vcf files
	VcfFile=$PairID".vcf"
	ConsequenceFile=$PairID"_csq.vcf"

	# Create a vcf with bcftools
	bcftools call $BcfFile --ploidy 1 --multiallelic-caller --output-type v --variants-only > $VcfFile

	# Progress
	echo -e "\e[1;31m Vcf file created, calling consequences... \e[0m"

	# Run bcftools csq to create a vcf file with the consequence calls
	bcftools csq --fasta-ref $PathToReference --gff-annot $PathToAnnotations $VcfFile -Ov -o $ConsequenceFile

	# Remove unecessary files
	rm $Samfile
	rm $Bamfile
	rm $SortedBam
	rm $IndexedBam
	rm $NDupBam
	rm $BcfFile
	rm $VcfFile

	# Progress
	echo -e "\e[1;31m Pair $PairID fully processed, moving onto next... \e[0m"

done
	
	
