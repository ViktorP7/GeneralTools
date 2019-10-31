# Load packages
library(ape)

# Set path variable for fasta file
fastaPath <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/sequences_Prox-10_26-09-2019.fasta"

# Read in the fasta file
sequences <- read.dna(fastaPath, as.character=TRUE, format = "fasta")

# Pull out sequences for animal one
animalOne <- sequences[c("16-5561_34.vcf.gz", "17-5559_50.vcf.gz"), ]
head(animalOne[, 1:10])

# Pull out sequences for animal two
animalTwo <- sequences[c("16-666_5.vcf.gz", "16-1297_4.vcf.gz", "16-1298_7.vcf.gz"), ]

# Function to check how many uniques present in a column
getUniqueAlleles <- function(fastaMatrix){
  
  # Get uniques for the given column
  uniques <- unique(fastaMatrix)
  
  # Get rid of the N positions from the count
  uniques <- uniques[!grepl("n",uniques)]
  
  # Get the length
  alleles <- length(uniques)
  
  # Return final
  return(alleles)
  
}


# Run the function using apply for the entire fasta for both animals
animalOneUniques <- apply(animalOne, 2, getUniqueAlleles)
animalTwoUniques <- apply(animalTwo, 2, getUniqueAlleles)

# Subset out any alleles which are greater than 1
animalOneAlleleIndexes <- which(animalOneUniques > 1)
animalTwoAlleleIndexes <- which(animalTwoUniques > 1)

# Plot a plot for animal one

# Plot a plot for animal two
plot(x=NA, y=NA, xaxt = "n", xlim = c(), ylim = c(0,2), main = "Troublesome Animal Two", 
     ylab = "SNP difference", las = 2)



