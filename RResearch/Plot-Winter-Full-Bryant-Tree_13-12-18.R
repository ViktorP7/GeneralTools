# 11/03/19 modified for a RAxML tree input file
# 25-03-19 modified for extra isolates

# Load packages
library(ape)
library(geiger)
library(phytools)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/CVRL-MAP-Batch-Jan.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_18-03-19"

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
plot.phylo(TheTree, edge.width = 0.2, font = 1, label.offset = 0.01, 
           show.tip.label = FALSE)
nodelabels()
tiplabels()

pdf("CombinedMarch2018Tree4.pdf", width=20, height=20)

# Root the tree at K10 
rootree <- root(TheTree, "MAPK10")

# Drop tips
dropNumbers <- c(61:74)
droppedTree <- drop.tip(rootree, dropNumbers)

# Drop the sheep
dropSheep <- c(64:77)
noSheepTree <- drop.tip(droppedTree, dropSheep)

# Drop the distant ones
dropDist <- c(60:63)
distDroppedTree <- drop.tip(noSheepTree, dropDist)

# Convert branch lengths to SNP values
distDroppedTree$edge.length <- distDroppedTree$edge.length * 49434

# Get the floored values o the lengths
flooredSNPs <- floor(distDroppedTree$edge.length)

# Get the colours
tipColours <- makeRegionColours(distDroppedTree$tip.label)

# Plot the no sheep tree
plot.phylo(distDroppedTree, edge.width = 0.2, font = 1, label.offset = 0.2, 
           tip.color = tipColours,
           align.tip.label = FALSE, type="phylogram", cex = 0.65)

# Add the SNP scale
add.scale.bar(cex = 3)

# Add SNP lengths to each edge
#edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), frame = "none", cex = 0.5)

dev.off()

# Function to get the labels names for bryant isolates
getBryantLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  chopVector <- TheTree$tip.label
  
  # Remove trailing _x.vcf.gz
  chopchopVector <- sapply(strsplit(chopVector, split = "_"), function(x) (x[1]))
  
  # Remove the > symbol
  nameVector <- sapply(strsplit(chopchopVector, split = ">"), function(x) (x[2]))
  
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
      }else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"],
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
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

# Function to generate colours based on species
makeSpeciesColours <- function(realNames){
  
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

# Function to generate colours based on region
makeRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if(grepl("UK", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Italy", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue" 
    } else if(grepl("Spain", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("France", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Scotland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("England", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Wales", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Germany", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Netherlands", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Czech", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Greece", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Norway", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("NewZealand", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("USA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Canada", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Venezuela", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("India", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Argentina", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("ERR0", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey"
    } else if(grepl("MAPK10", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey"
    } else {
      
      colourVec[index] <- "green"
    }
  }
  
  return(colourVec)
}
# To get the scale bar, times all the edge lengths by the length of sequences
# This can be found in the fasta file at the top
# Example - rooted$edge.length <- rooted$edge.length * 671
# Then use the add.scale.bar() command
