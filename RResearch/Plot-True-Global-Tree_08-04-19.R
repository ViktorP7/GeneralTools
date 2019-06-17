# 11/03/19 modified for a RAxML tree input file
# 25-03-19 modified for extra isolates
# 08-04-19 modified for better root and to find VNTR group distances
# 12-06-19 modified for new influx of data

# Load packages
library(ape)
library(phytools)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MAP-Metadata-Formatted-May19.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_10-06-19"

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
           show.tip.label = TRUE, cex = 0.1)
nodelabels(cex = 0.05, frame = "none")
tiplabels()

pdf("CombinedJune2019Trees.pdf", width=20, height=20)

# Root the tree at 466 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 466)

# Drop tips for far away ancestors (silvaticum and hominissius)
dropNumbers <- c(186,187)
droppedTree <- drop.tip(rootree, dropNumbers)

# Extract the clade that doesn't have all the distant sheep and cows
extractedTree <- extract.clade(droppedTree, node = 452)

# Convert branch lengths to SNP values
extractedTree$edge.length <- extractedTree$edge.length * 48283

# Get the rounded values o the lengths
roundedSNPs <- round(extractedTree$edge.length)

# Assign rounded SNPs
extractedTree$edge.length <- roundedSNPs

# Get the colours
tipColours <- makeRegionColours(extractedTree$tip.label)

# Plot the tree
plot.phylo(extractedTree, edge.width = 0.2, font = 1, label.offset = 0.2, 
           tip.color = tipColours,
           align.tip.label = FALSE, type="phylogram", cex = 0.2)

# Add the SNP scale
add.scale.bar(cex = 3)

# Plot international tree as a fan
plot.phylo(extractedTree, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = tipColours, edge.color = "grey50",
           align.tip.label = FALSE, type="fan", cex = 0.5, show.tip.label = FALSE)

#Add shaped tip labels
tiplabels(pch = 18, col = tipColours,  cex = 2)

# Add the SNP scale
add.scale.bar(x=30,y=-120, cex = 3,lcol = "grey50", lwd = 3)
text(x=60,y=-130, "SNPs", cex = 3)

# Add a legend
legend(x=40, y=140, legend = c("Ireland", "Europe", "World"), 
       text.col = c("darkgreen", "blue", "darkorange3"), bty = "n", cex = 2,y.intersp = 0.6)

dev.off()

#### Functions ####

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
        
        nameVector[index] <- "14-6278_Cork_9_1"
        
      }else if(nameVector[index] == "14-5154"){
        
        nameVector[index] <- "16-5154"
        
      }else if(nameVector[index] == "16-4434"){
        
        nameVector[index] <- "16-4934"
        
      }else if(nameVector[index] == "14-4776"){
        
        nameVector[index] <- "17-4776"
        
      }else if(nameVector[index] == "14-2662"){
        
        nameVector[index] <- "14-2622"
      }else if(nameVector[index] == "14-7468"){
        
        nameVector[index] <- "14-7486"
        
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
      if(isoTable[row,"AliquotFormat"] == nameVector[index]){
        
        newname <- paste(nameVector[index], "_", isoTable[row,"Herd Ref"], "_",
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

# Function to generate colours based on region
makeRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkgreen"
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
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("USA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("Canada", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("Venezuela", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("India", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("Argentina", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange3"
    } else if(grepl("ERR0", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey30"
    } else if(grepl("MAPK10", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey30"
    } else {
      
      colourVec[index] <- "darkgreen"
    }
  }
  
  return(colourVec)
}