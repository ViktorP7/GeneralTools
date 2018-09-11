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
nodelabels(cex = 0.4)

pdf("tree.pdf", width=20, height=20)

# Root the tree at node 205 - like in the paper
tree <- root(TheTree, node=205)
nodelabels(cex = 0.4)

# Rotate the clades 
rotTree <- rotate(tree, node = 148, polytom = c(2,3))

# Plot overview tree
plot.phylo(rotTree, edge.width = 0.2, font = 1, label.offset = 0.001,
           align.tip.label = TRUE, type="phylogram", show.tip.label=FALSE)

# Remove tips for M. avium and silvaticum isolates, and other ERR's
dropNumbers <- c(57,58,62,115,80,81)
droppedTree <- drop.tip(rotTree, dropNumbers)

# Get the colours
tipColours <- makeColours(droppedTree$tip.label)

# Convert branch lengths to SNP values
droppedTree$edge.length <- droppedTree$edge.length * 39423

# Get the floored values o the lengths
flooredSNPs <- floor(droppedTree$edge.length)

# Plot the new tree
plot.phylo(droppedTree, edge.width = 0.2, font = 1, label.offset = 0.001,
           cex=0.7, align.tip.label = TRUE, type="phylogram", tip.color = tipColours)
legend("left", legend=c("Colouring based on host","Species.Location.IsolateNo.INMV", 
                        "Human", "Deer", "Cattle/Bovidae", "Sheep/Capridae",
                        "Passaged", "Other", "Scale shows no. SNPs"), text.col = c("Black","Black", "Red",
                                                           "Orange", "Green",
                                                           "Blue", "Grey", "Black", "Black"),
       bty = "n", cex = 1)

add.scale.bar()

# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), bg = "yellow", cex = 0.4)

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
      
      colourVec[index] <- "red"
    } else if(grepl("Cow", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
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

# Position scale line between last two major x-axis tick marks
# and 1/10th of the total y-range above the lower y-axis coordinate
#lines(c(floor(axisLimits[2]-1),floor(axisLimits[2])),     
#      rep(axisLimits[3] + 0.1*(axisLimits[4] - axisLimits[3]), 2))

# Place the units label at the midpoint of and just below the scale line
#text(x=mean(c(floor(axisLimits[2]), floor(axisLimits[2]))), 
#     y=axisLimits[3] + 0.1*(axisLimits[4] - axisLimits[3]),
#     label="SNPs", adj=c(0.5, 1.5))