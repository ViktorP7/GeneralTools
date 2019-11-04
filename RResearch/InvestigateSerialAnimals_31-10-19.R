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
plot(x=NA, y=NA, xaxt = "n", xlim = c(as.Date("2016-09-01"), as.Date("2017-11-01")), ylim = c(0,6), main = "Troublesome Animal One", 
     ylab = "SNP difference", xlab = "Time", las = 2)
axis.Date(side = 1, x=c("2016-09-01", "2017-11-01"), format = "%Y-%m")
lines(x=c(as.Date("2016-09-01"),as.Date("2017-11-01")) , y = c(0,0), col = "purple")
lines(x=c(as.Date("2017-09-04"),as.Date("2017-09-04")) , y = c(0,5), col = "red")
lines(x=c(as.Date("2017-09-04"),as.Date("2017-11-01")) , y = c(5,5), col = "red")
points(x=c(as.Date("2016-10-07"),as.Date("2017-09-04")), y=c(0,5), 
       col=c("blue","red"), cex = 2, pch = 19)
text(x=c(as.Date("2016-10-07"),as.Date("2017-09-04")), y = c(0,5),
     labels = c("16-5561", "17-5559"), cex = 0.7, pos = 3)

# Plot a plot for animal two
plot(x=NA, y=NA, xaxt = "n", xlim = c(as.Date("2016-02-01"), as.Date("2016-03-10")), ylim = c(0,3), main = "Troublesome Animal Two", 
     ylab = "SNP difference", xlab = "Time", las = 2)
axis.Date(side = 1, x=c("2016-02-01", "2016-03-10"), format = "%Y-%m-%d")
lines(x=c(as.Date("2016-02-01"),as.Date("2016-03-10")) , y = c(0,0), col = "purple")
lines(x=c(as.Date("2016-03-03"),as.Date("2016-03-03")) , y = c(0,2), col = "red")
lines(x=c(as.Date("2016-03-03"),as.Date("2016-03-10")) , y = c(2,2), col = "red")
points(x=c(as.Date("2016-02-02"),as.Date("2016-03-03"),as.Date("2016-03-03")), y=c(0,0,2), 
       col=c("blue","purple","red"), cex = 2, pch = 19)
text(x=c(as.Date("2016-02-02"),as.Date("2016-03-03"),as.Date("2016-03-03")), y = c(0,0,2),
     labels = c("16-666", "16-1298 Faeces", "16-1297 Ileocaecal"), cex = 0.7, pos = 3)




