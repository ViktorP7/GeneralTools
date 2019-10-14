
# Sequence Processing Pipeline and Visualisation 14-10-19


## Installing necessary programs to Virtual Box
Hopefully you've had no trouble installing R, R Studio and Virtual Box. The next thing to do is to create an Ubuntu virtual machine in Virtual Box and set up a shared folder. The instructions for doing this are available here: https://josephcrispell.github.io/BlogPosts/VirtualBoxUbuntu_09-11-17/VirtualBoxUbuntu_09-11-17.html

Once this is set up, we need to install several tools, some of which you've downloaded already. The instructions for installing these tools are available here: https://github.com/JosephCrispell/GeneralTools/blob/master/ProcessingPipeline/installation.txt

To properly install some of these tools, it may be necessary to install a couple of dependencies first. Type the following code in the Ubuntu terminal:

```bash

# Install dependencies
sudo apt-get install libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev

# Install Java, pigz and raxml
sudo apt install default-jre pigz raxml
```

Finally, go to: https://github.com/JosephCrispell/GeneralTools and download the repository, and unzip it into your shared folder. You should now have everything in your shared folder to run the bioinformatics pipeline.

## Running the bioinformatics pipeline
The general commands for running the pipeline can be found here: https://github.com/JosephCrispell/GeneralTools/blob/master/ProcessingPipeline/README.txt 

However, we'll run the specific commands printed below to make things easier. First, navigate to your folder with all the fastq.gz format sequence files.

```bash
# Change directory to sequence files folder
cd SequenceFiles
```

Check what is present in the folder by typing <code>ls</code>. The first thing we need to do is run a quality check on all the sequence files using FastQC, which will enable us to determine trimming parameters. 

```bash
# Run fastqc on all fastq.gz files
fastqc *fastq.gz
```

Instead of running a piece of code again and again for each of the paired files, we can instead use <code>*fastq.gz</code> to denote all files ending in fastq.gz, saving us a lot of time. You will see the asterisk symbol being used again to denote all files of a given file ending. 

Once FastQC is finished, our folder will look like a mess, so we should tidy it up a bit. Make a folder for the HTML files. We then want to move all our HTML files into this folder.

```bash
# Make a folder called HTMLs
mkdir HTMLs

# Move all .html files into the folder
mv *.html HTMLs
```

This leave us with a lot of .zip files, which we do not need. We can remove them using <code>rm</code>. Be very careful using this command as it can delete EVERYTHING on an Ubuntu machine if used improperly!

```bash
# Remove all zip files
rm *.zip
```

Now, let's have a look at all the .html files for all the paired files. You should make a log file and note down the characteristics of each pair. They should all generally conform to the same quality pattern, but if there are any oddities such as low %GC (you're expecting about 68% for MAP), note it down! The template CutAdapt settings file should already contain appropriate trimming parameters, however if you think any need changing, go ahead and do so. 

Once you're satisfied, it's time to run the ProcessRawReads script. This is a bash script that acts as a pipeline. A pipeline runs files throught one program, then takes the outputs and runs them through another and so on. The benefit of having a pipeline is that you can do a complicated series of steps by typing just one command into your terminal!

In this case we are passing our fastq.gz files into CutAdapt to trim the reads as per our settings file. The trimmed reads are then passed on to bwa mem, which aligns the reads to the reference genome. Any reads that do not align to the reference are put through BLAST to determine their origin. The sequence alignment map (SAM) is converted between a couple of different formats via samtools, eventually being converted into a variant calling file (VCF) with bcftools. The VCF format shows the entire genome, whether an alternate allele is present, and all the relevant quality info.

Next, move into the vcf containing folder, where we want to merge the VCF files into one files, making it easier for subsequent scripts to work.

```bash
# Process raw reads
bash ../JoeGeneralTools/ProcessingPipeline/ProcessRawReads_28-06-17.sh fastq.gz cutadapt ../Reference/MAPK10.fasta ../JoeGeneralTools/ProcessingPipeline/PickRandomReadsFromSAM_13-07-17.pl ../JoeGeneralTools/ProcessingPipeline/ExamineBLASTOutput_13-07-17.pl

# Move to vcf folder
cd vcfFiles

# Run merge tool
java -jar ../../JoeGeneralTools/ProcessingPipeline/MergeVCFFiles_10-01-19.jar . ../../Reference/MAPK10.gff3
```

We can check the genome coverage of the isolates. After creating the merged VCF, we need to filter it to get rid of bad quality or badly supported positions.

```bash
# Check genome coverage
perl ../../JoeGeneralTools/ProcessingPipeline/ExamineGenomeCoverage_03-07-17.pl 20 genomeCoverage_16-10-19.txt

# Filter the merged file
perl ../../JoeGeneralTools/ProcessingPipeline/FilterVariants_15-09-17.pl 1 30 4 35 0.95 0.05 0 0 merged_16-10-19.txt
```

Now we want to check the isolate coverage. Add the prefix "PreRescue" to the name of the output isolate coverage file. If some isolates are looking poor, we can run the rescue script after, which will look at every variant position and attempt to "rescue" depending on whether or not it's present in other isolates.

```bash
# Calculate isolate coverage
perl ../../JoeGeneralTools/ProcessingPipeline/CalculateIsolateCoverageFromFiltered_20-09-17.pl filtered_30-4-35-0.95-0.05-0-0_16-10-19.txt

# Run rescue tool
perl ../../JoeGeneralTools/ProcessingPipeline/RescueVariantPositionInfo_18-09-17.pl 3 5 0.95 filtered_30-4-35-0.95-0.05-0-0_16-10-19.txt
```

Now, we can examine the isolate coverage on the rescued file and see if the rescue has helped any of our poorer quality isolates. Make sure to name the output file "PostRescue".

```bash
perl ../../JoeGeneralTools/ProcessingPipeline/CalculateIsolateCoverageFromFiltered_20-09-17.pl filtered-rescued_3-5-0.95_16-10-19.txt
```

Finally we need to create a FASTA file from our filtered variants. There is a proximity filter used here to get rid of variants within a given short distance of each other, as they are probably the same thing.

```bash
perl ../../JoeGeneralTools/ProcessingPipeline/CreateFastaWithReferenceFromFiltered_28-06-17.pl 1 10 ../../Reference/MAPK10.fasta filtered-rescued_3-5-0.95_16-10-19.txt
```

## Treebuilding
We use a program called RAxML for treebuilding. It takes our FASTA file as input and uses the maximum likelihood algorithm to construct the tree. We can also designate a bootstrap replicate value. While it can be launched via command line, we will use the following R script to launch RAxML instead, with some changes: https://github.com/JosephCrispell/GeneralTools/blob/master/RunRAXML_27-06-18.R

```r
# Note the FASTA file name and path
fastaFile <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SequenceFiles/vcfFiles/sequences_Prox-10_14-10-2019.fasta"

# Set the working directory
setwd("C:/Users/UCD/Desktop/UbuntuSharedFolder/SequenceFiles/")

# Build analysis name
analysisName <- "RaxML-R_14-10-19"

# Set the input parameters for RAXML
model <- "GTRCAT" # No rate heterogenity
seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
nBootstraps <- 100
nThreads <- 2

# Build the command for RAXML
command <- paste("raxmlHPC-PTHREADS", # Note on WINDOWS replace with: /path/to/raxmlHPC-PTHREADS.exe
                 " -f a", # Algorithm: Rapid boostrap inference
                 " -N ", nBootstraps,
                 " -T ", nThreads,
                 " -m ", model, " -V", # -V means no rate heterogenity
                 " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                 " -n ", analysisName,
                 " -s ", fastaFile, sep="")

# Run RAXML
system(command, intern=TRUE)
```

## Tree Visualisation in R
Now that we have made our tree file, we need to import it into R along with the associated metadata file, which will contain matching isolate IDs and VNTR info for each. Before we do that, we need to install a number of R packages.

```r
# Install relevant packages
install.packages("ape", "phytools", "scales")

# Load packages
library(ape)
library(phytools)
library(scales)

# Define paths to tree and metadata
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SequenceFiles/RAxML_bipartitions.RaxML-R_16-10-19"
pathData <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SequenceFiles/Chuhan_Metadata.csv"

# Read in metadata
metadataTable <- read.table(pathData, header = TRUE, sep = ",",
                            stringsAsFactors=FALSE,
                            check.names=FALSE)
                            
# Read in tree
theTree <- read.tree(pathTree)
```

We've now got the table and tree read into R, but now we need to combine the two sets of information somehow. In the tree file, if we interrogate <code>theTree$tip.label</code>, you'll notice that the tip labels of the tree feature the isolate numbers. We need to match these to their metadata entries and attach on VNTR info to the tip labels. We can do this by writing a custom function. We've used a couple of functions already to read in the tree and table. Think of a function as a machine - you give it input(s), it does something to them and returns an output. Let's write a function to attach on VNTR info to our tip labels.

```r
# Name the function and define inputs
addVNTRtoTipLabel <- function(tree, data){

  # Create a vector to match the tip label vector in the tree
  nameVector <- tree$tip.label
  
  # Loop thru the data
  for(row in 1:nrow(data)){
  
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
    
      # Check if the current name matches the tip label, and if so...
      if(data[row, "Name"] == nameVector[index]){
      
        # Paste together the label with the VNTR number
        newLabel <- paste(nameVector[index], "_", data[row, "VNTR"])
        
        # Store the label in the corresponding place in the name vector
        nameVector[index] <- newLabel
      
      # But if the names don't match, then...
      }else{ 

        # Skip to next iteration
        next
        
      }
    }
  }
 
 # Return an output from the function
 return(nameVector)
}
```

There are a lot of things happening in this function, don't worry if you don't understand the specifics. The key thing to follow is that the function works by matching tip label names to names in the metadata file, and then returns the tip labels with the VNTR number attached on to the end of each label. Now that the function is written, it needs to be run so that it's loaded into memory, and then we can use it to get the tip labels we want and place them in the tree.


