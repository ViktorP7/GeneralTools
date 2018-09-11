# Load packages
library(ape)
library(geiger)
library(phytools)

# Set path variables
pathIso <- "path/to/csv"
path <- "path/to/tree"

# Read in table of isolates
isoTable <- read.table(pathIso,
                       header = TRUE,
                       sep = ",",
                       stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted

# Read in tree
TheTree <- read.tree(path)

# Get the new names
realNames <- getLabels(isoTable, TheTree)

# Change the names in the tree
TheTree$tip.label <- realNames

# Plot the tree
plot.phylo(TheTree, edge.width = 0.2, font = 1, label.offset = 0.01)
nodelabels()
tiplabels()

pdf("irishtree.pdf", width=20, height=20)

# Root the tree at node 205 - like in the paper
rootree <- root(TheTree, "MAPK10")

# Get the colours
tipColours <- makeColours(rootree$tip.label)

# Convert branch lengths to SNP values
rootree$edge.length <- rootree$edge.length * 671

# Get the floored values o the lengths
flooredSNPs <- floor(rootree$edge.length)

# Plot the new tree
plot.phylo(rootree, edge.width = 0.2, font = 1, label.offset = 0.001,
          align.tip.label = TRUE, type="phylogram", tip.color = tipColours, cex = 1.5)
legend("bottomleft", legend=c("Species.Location.IsolateNo.INMV"), text.col = c("Black"),
       bty = "n", cex = 1.5)

# Add the SNP scale
add.scale.bar(cex = 1.5)

# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), bg = "yellow", cex = 1.5)
edgelabels("SNPs", adj = c(0.5, 1.25), bg = "lightblue", cex = 1.5)

dev.off()

# Function to get the labels names
getLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"], 
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
        nameVector[index] <- newname
      } else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"], 
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
        nameVector[index] <- newname
      } else if(nameVector[index] == "Ref-1997") {
        
        nameVector[index] <- "MAP K10"
      } else{
        
        next
      }
      
    }
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to generate colours based on species
makeColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Human", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Cow", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Bison", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green" 
    } else if(grepl("Sheep", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Goat", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Moufflon", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Deer", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("Passaged", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey"
    } else {
      
      colourVec[index] <- "black"
    }
  }
  
  return(colourVec)
}

# To get the scale bar, times all the edge lengths by the length of sequences
# This can be found in the fasta file at the top
# Example - rooted$edge.length <- rooted$edge.length * 671
# Then use the add.scale.bar() command
