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
plot.phylo(TheTree, edge.width = 0.2, font = 1, label.offset = 0.01, cex = 0.4)
#nodelabels(cex = 0.4)
#tiplabels(cex = 0.4)

droppedTree <- drop.tip(TheTree, 142)
# Plot the tree
plot.phylo(droppedTree, edge.width = 0.2, font = 1, label.offset = 0.01, cex = 0.4)

pdf("ahlstromtree.pdf", width=20, height=20)

# Root the tree at node 205 - like in the paper
rootree <- root(droppedTree, "MAPK10")

# Convert branch lengths to SNP values
rootree$edge.length <- rootree$edge.length * 29366

# Get the floored values o the lengths
flooredSNPs <- floor(rootree$edge.length)

# Plot the new tree
plot.phylo(rootree, edge.width = 0.2, font = 1, label.offset = 0.001,
           align.tip.label = TRUE, type="phylogram")


# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), bg = "yellow", cex = 0.5)

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
        
        nameVector[index] <- isoTable[row, "UCID"]
      
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
