# 24-04-19 Reworked script for updated inputs and better tree

# Load packages
library(ape)
library(phytools)
library(scales)

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
           show.tip.label = TRUE, cex = 0.1)
nodelabels(cex = 0.05, frame = "none")
tiplabels()

# Root the tree at 332 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 332)

# Drop tips for far away ancestors (silvaticum and hominissius)
dropNumbers <- c(69,70)
droppedTree <- drop.tip(rootree, dropNumbers)

# Extract the clade that doesn't have all the distant sheep
extractedTree <- extract.clade(droppedTree, node = 319)

# Extract further
extractTree <- extract.clade(extractedTree, node = 227)

# Drop the ERR0s
dropErr <- c(57,58,100)
noErrTree <- drop.tip(extractTree, dropErr)

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(noErrTree$tip.label)
irishOnlytree <- drop.tip(noErrTree, dropper)

# Extract CIT region
CITTree <- extract.clade(irishOnlytree, node = 134)

# Convert branch lengths to SNP values
CITTree$edge.length <- CITTree$edge.length * 49434

# Get the floored values o the lengths
flooredSNPs <- floor(CITTree$edge.length)

# Get the colours
tipColours <- makeIrishRegionColours(CITTree$tip.label)

# Plot the national tree
plot.phylo(CITTree, edge.width = 0.2, font = 1, label.offset = 0.2, 
           tip.color = tipColours,
           align.tip.label = FALSE, type="phylogram", cex = 0.7)


# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), 
           bg = "white", cex =1, frame = "n")


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
        
        nameVector[index] <- "14-6278_Cork 9(Jul18)_1"
        
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

# Function to generate colours based on region
makeIrishRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("NIreland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Louth", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if(grepl("Cavan", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red" 
    } else if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("Leitrim", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Clare", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("Wexford", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if(grepl("Limerick", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("Wicklow", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if(grepl("Cork", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("NA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Kerry", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if(grepl("Donegal", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("CIT", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if (grepl("MAPK10", colourVec[index]) == TRUE){
      
      colourVec[index] <- "purple"
    } else if (grepl("Dublin", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if (grepl("Laois", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if (grepl("Kilkenny", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if (grepl("Meath", colourVec[index]) == TRUE){
      
      colourVec[index] <- "green"
    } else if (grepl("Galway", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if (grepl("Tipperary", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orange"
    } else if (grepl("Roscommon", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else {
      
      colourVec[index] <- "black"
    }
    
  }
  
  return(colourVec)
}

# Function to create vector with international tips to drop
toDropInternationalTips <- function(tiplabel){
  
  # Create vector to store index values of tips to be dropped
  dropVector <- c()
  
  # Loop thru tip labels and drop as required
  for(index in 1:length(tiplabel)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("UK", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Italy", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index) 
    } else if(grepl("Spain", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("France", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Scotland", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("England", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Wales", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Germany", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Netherlands", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Czech", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Greece", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Norway", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("NewZealand", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("USA", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Canada", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Venezuela", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("India", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Argentina", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("ERR0", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } 
  }
  return(dropVector)
}
