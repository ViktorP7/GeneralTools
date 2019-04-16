# Load packages
library(ape)
library(geiger)
library(phytools)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/CVRL-MAP-Batch.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/WinterMAPTree.tree"

# Read in table of bryant isolates
isoBryantTable <- read.table(pathBryantIso,
                       header = TRUE,
                       sep = ",",
                       stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted

# Read in table of CVRL isolates
isoCVRLTable <- read.table(pathNewIso,
                             header = TRUE,
                             sep = ",",
                             stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                             check.names=FALSE) # Names left as they are, no dots inserted

# Read in tree
TheTree <- read.tree(pathTree)

# Get the Bryant names
realNames <- getBryantLabels(isoBryantTable, TheTree)

# Update the names in the tree
TheTree$tip.label <- realNames

# Get metadata for CVRL isolates
realNames <- getCVRLLabels(isoCVRLTable, TheTree)

# Update the names in the tree
TheTree$tip.label <- realNames

# Plot the tree
plot.phylo(TheTree, edge.width = 0.2, font = 1, label.offset = 0.01)
nodelabels()
tiplabels()

pdf("December2018Tree.pdf", width=20, height=20)

# Root the tree at K10 
rootree <- root(TheTree, "MAPK10")

# Convert branch lengths to SNP values
rootree$edge.length <- rootree$edge.length * 1789

# Get the floored values o the lengths
flooredSNPs <- floor(rootree$edge.length)

# Plot the new tree
plot.phylo(rootree, edge.width = 0.2, font = 1, label.offset = 0.001,
           align.tip.label = TRUE, type="phylogram", cex = 0.9)

# Add the SNP scale
add.scale.bar(cex = 1.0)

# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), frame = "none", cex = 0.5)

dev.off()

# Function to get the labels names for bryant isolates
getBryantLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"], "(Bryant)")
        nameVector[index] <- newname
      }else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"], "(Bryant)")
        nameVector[index] <- newname
      }else if(nameVector[index] == "Ref-1997") {
        
        nameVector[index] <- "MAP K10"
      }else if(nameVector[index] == "6287-MAP"){
        
        nameVector[index] <- "14-6278_Cork 9_1(July2018)"
        
      }else if(nameVector[index] == "14-5154"){
        
        nameVector[index] <- "16-5154"
        
      }else if(nameVector[index] == "16-4434"){
        
        nameVector[index] <- "16-4934"
        
      }else{
        
        next
      }
      
    }
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to get the labels names fr CVRL isolates
getCVRLLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Isolate-ID(TBxx-xxxxxx)"] == nameVector[index]){
        
        newname <- paste(nameVector[index], "_", isoTable[row,"Herd Location"], "_",
                         isoTable[row,"INMV Group"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}


# To get the scale bar, times all the edge lengths by the length of sequences
# This can be found in the fasta file at the top
# Example - rooted$edge.length <- rooted$edge.length * 671
# Then use the add.scale.bar() command
