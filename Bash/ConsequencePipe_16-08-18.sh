#!/bin/bash

# Trimming, alignment, and consequence calling of paired raw reads
# Contains adapted code from ProcessRawReads_09-08-18 by Joe Crispell
# Author - Viktor Perets 16/08/18
# For bcftools csq to work properly, your fasta and annotations should be located in the working directory

# Command structure
# ConsequencePipe_16-08-18.sh FastqFileEnding PathToPrinseq PathToReference PathToAnnotations

# Define fastq file ending variable, path to prinseq, reference and annotations
FastqFileEnding=$1
PathToPrinseq=$2
PathToReference=$3
PathToAnnotations=$4

# Define prinseq settings
LENGTH=75			# The minimum length of read to be accepted
MEANQUAL=25			# Filter sequence if mean quality score below x
TRIML=15			# Trim sequence at the 5' end by x positions
TRIMR=5				# Trim sequence at the 3' end by x positions
TRIMQUALL=25		# Trim sequence by quality score from the 5' end with this threshold
TRIMQUALR=25		# Trim sequence by quality score from the 3' end with this threshold score
TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
WINDSIZE=10			# The sliding window size used to calculate quality score by type
TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
TRIMRTAIL=5			# Trim poly A/T > X length at 3' end

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
	echo -e "\e[1;31m Unzipping raw reads of $PairID... \e[0m"
	
	# Unzip each file
	gunzip $PairOne
	gunzip $PairTwo

	# Remove .gz from file names
	PairOne=`echo ${PairOne:0:-3}`
	PairTwo=`echo ${PairTwo:0:-3}`

	# Progress
	echo -e "\e[1;31m Files unzipped, beginning trimming... \e[0m"

	# Trim the reads using prinseq
	perl $PathToPrinseq -fastq $PairOne -fastq2 $PairTwo -min_len $LENGTH -min_qual_mean $MEANQUAL -trim_left $TRIML -trim_right $TRIMR -trim_qual_left $TRIMQUALL -trim_qual_right $TRIMQUALR -trim_qual_type $TRIMTYPE -trim_qual_window $WINDSIZE

	# Progress
	echo -e "\e[1;31m Trimming complete, zipping the files back up... \e[0m"
	
	# Zip the files back up
	pigz $PairOne
	pigz $PairTwo

	# Progress
	echo -e "\e[1;31m Files zipped up, beginning alignment of trimmed files... \e[0m"

	# Find output files and store them 
	TrimmedFiles=(`ls | grep $PAIRID".*prinseq_good" | grep -v "singletons"`)
	PairOne=${TrimmedFiles[0]}
	PairTwo=${TrimmedFiles[1]}
	
	# Remove unecessary files
	rm *prinseq_bad*
	rm *prinseq_good_singletons*

	# Variables to store sam, bam sorted bam and indexed bam
	Samfile=$PairID"_aln-pe.sam"
	Bamfile=$PairID".bam"
	SortedBam=$PairID"_srtd.bam"
	IndexedBam=$SortedBam".bai"
	NDupBam=$PairID"_srtd_ndup.bam"

	# Align to reference with bwa mem and pipe to samtools to create, sort and index a bam file
	bwa mem $PathToReference $PairOne $PairTwo > $Samfile
	
	# Remove Files
	rm $PairOne
	rm $PairTwo
	
	samtools view -b $Samfile > $Bamfile
	samtools sort $Bamfile -o $SortedBam
	samtools index $SortedBam
	samtools rmdup $SortedBam $NDupBam

	# Progress
	echo -e "\e[1;31m Alignment complete, creating bcf file... \e[0m"

	# Variable to store bcf file
	BcfFile=$PairID".bcf"

	# Create a bcf file with samtools
	samtools mpileup --adjust-MQ 5 --min-MQ 5 --min-BQ 5 --uncompressed --fasta-ref $PathToReference $NDupBam > $BcfFile

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

	# Zip the vcf files up
	pigz $VcfFile

	# Progress
	echo -e "\e[1;31m Pair $PairID fully processed, moving onto next... \e[0m"

done
	
	
