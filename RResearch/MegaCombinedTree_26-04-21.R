# Created 18-02-20
# Script to process and visualise phylogenetic trees
# Please run this FIRST, as it is a PREREQUISITE for the other scripts to work

# Load packages
library(ape)
library(phytools)
library(scales)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaOct2020Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathBigTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/GenieBigMAP/RAxML_bipartitions.RaxML-R_23-04-21"
pathNI <- "C:/Users/UCD/Documents/Lab/CVRL MAP/NIMetaOct2020.csv"
pathAus <- "C:/Users/UCD/Documents/Lab/AustraliaMetadata.csv"
pathCan <- "C:/Users/UCD/Documents/Lab/CanadaMetadata.csv"

# Read in table of bryant isolates
isoBryantTable <- read.table(pathBryantIso, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of CVRL isolates
isoCVRLTable <- read.table(pathNewIso, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of NI isolates
isoNITable <- read.table(pathNI, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of Aus isolates
ausTable <- read.table(pathAus, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

canTable <- read.table(pathCan, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

shortCounties <- c("AnNI", "ArNI", "CW", "CN", "CE", "C", "DL", "DoNI", "D", "FNI", "G", "KY", "KE", "KK", "LS", "LM",
                   "L", "DeNI", "LD", "LH", "MO", "MH", "MN", "OY", "RN", "SO", "T", "TNI", "W", "WH", "WX", "WW")

#### Tree file processing ####

# Read in tree
BigTree <- read.tree(pathBigTree)

# Root on another node
BigRoot <- root(BigTree, node = 1292)

# Drop sheep isolates and avium
tipstogo <- c(457:556)
noSTree <- drop.tip(BigRoot, tipstogo)

# Get the Bryant names
realNames <- newgetBryantLabels(isoBryantTable, noSTree)

# Update the names in the tree
noSTree$tip.label <- realNames

# Get metadata for CVRL isolates
realNames <- getCVRLLabels(isoCVRLTable, noSTree)

# Update the names in the tree
noSTree$tip.label <- realNames

# Get metadata for NI isolates
realNames <- getNILabels(isoNITable, noSTree)

# Update names
noSTree$tip.label <- realNames

# GET CAN metadata
realNames <- getSRALabels(canTable, noSTree)

noSTree$tip.label <- realNames

# Get Aus metadata
realNames <- getSRALabels(ausTable, noSTree)

noSTree$tip.label <- realNames

# Convert branch lengths to SNP values
noSTree$edge.length <- round(noSTree$edge.length * 31255)

worldDist <- cophenetic(noSTree)

# Make EU colours
worldCols <- newmakeRegionColours(noSTree$tip.label)

# Make EU plot
outputFile <- paste("EU-Tree_29-01-21.png", sep="")
png(outputFile, height=5000, width=4500)

# Plot EU tree
plot.phylo(euOnlyTree, edge.width = 10, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(505,456,418,402,386,373,351,305), frame = "n", cex=15, adj = c(1,0), col = "red")



#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 10)

# Add the SNP scale
add.scale.bar(x=0, y=0, cex = 8, lwd = 15)
text(x=65, y=0, cex = 8, "SNPs")

# Add a legend
legend(x=150, y=300, legend = c("Ireland", "UK", "England", "Scotland", "Wales",
                                "Italy", "Spain", "France", "Germany", "Netherlands",
                                "Czech Rep.", "Greece", "Norway"), 
       text.col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
                    "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
                    "mediumblue", "slateblue", "purple"), pch = 18,
       col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
               "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
               "mediumblue", "slateblue", "purple"),
       bty = "n", cex = 8.8, y.intersp = 0.8, title = "Country")

dev.off()


#### Tree plotting (.pdf) ####


# Make EU plot pdf
outputFile <- paste("World-Tree_26-04-21.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot EU tree
plot.phylo(noSTree, edge.width = 6, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)


#Add shaped tip labels
tiplabels(pch = 18, col = worldCols,  cex = 6)

# Add the SNP scale
add.scale.bar(x=75, y=0, cex = 10, lwd = 15)
text(x=142, y=0, cex = 10, "SNPs")

# Add a legend
#legend(x=150, y=280, legend = c("Ireland", "UK", "England", "Scotland", "Wales",
#                                "Italy", "Spain", "France", "Germany", "Netherlands",
#                                "Czech Rep.", "Greece", "Norway"), 
#       text.col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
#                    "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
#                    "mediumblue", "slateblue", "purple"), pch = 18,
#       col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
#               "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
#               "mediumblue", "slateblue", "purple"),
#       bty = "n", cex = 10, y.intersp = 0.8, title = "Country", title.col = "black")

dev.off()
#### Functions ####

# Function to get the labels names for bryant isolates
newgetBryantLabels <- function(isoTable, TheTree){
  
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
      }else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row, "Host"],
                         isoTable[row,"Country of origin"], "_", splitter, "_",
                         isoTable[row,"INMV"])
        nameVector[index] <- newname
      }else if(nameVector[index] == "Ref-1997") {
        
        nameVector[index] <- "MAP K10"
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
        
        herd <- strsplit(isoTable[row,"Herd Identifier"], split = " ")[[1]][2]
        
        newname <- paste(nameVector[index], "_", isoTable[row, "County"], "_", herd, "_",
                         isoTable[row,"INMV Group"], "_", isoTable[row, "Birth County"], "_", isoTable[row, "Number of moves to herd of sampling"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to get labels off an SRA document
getSRALabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Run"] == nameVector[index]){
        
        newname <- paste(isoTable[row, "Collection_Date"], "_", isoTable[row, "geo_loc_name"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
}

# Function to get labels for NI isolates
getNILabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"SeqRef"] == nameVector[index]){
        
        county <- strsplit(isoTable[row,"Herd Ref"], split = " ")[[1]][1]
        herd <- strsplit(isoTable[row,"Herd Ref"], split = " ")[[1]][2]
        
        newname <- paste(isoTable[row,"AliquotFormat"], "_", county, "_", herd, "_",
                         isoTable[row,"INMV"], "_", isoTable[row, "Birth County"], "_", isoTable[row, "Herd of Birth"])
        nameVector[index] <- newname
      }else{
        
        next
      }
      
    }
    
  }
  nameVector <- gsub(" ", "", nameVector, fixed = TRUE)
  
  return(nameVector)
  
}

# Function to pull out matrix names and simplify to chosen
getNames <- function(mat, chosen){
  
  # Create a vector for names
  rowcolNames <- rep(NA, length(colnames(mat)))
  
  if(chosen == "VNTR"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(rownames(mat)[index], split = "_")[[1]][4]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  } else if(chosen == "Herd"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      one <- strsplit(colnames(mat)[index], split = "_")[[1]][2]
      two <- strsplit(colnames(mat)[index], split = "_")[[1]][3]
      
      vRow <- paste(one,"_",two)
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
    
  } else if(chosen == "CCounty"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][2]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
    
  } else if(chosen == "BCounty"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][5]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  } else if(chosen == "Same"){
    
    # Loop thru and split off what's needed and fill into vectors
    for(index in 1:length(rowcolNames)){
      
      # Store row name
      vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][6]
      
      rowcolNames[index] <- vRow
    }
    return(rowcolNames)
  }
}

# Function to generate colours based on region
newmakeRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkgreen"
    } else if(grepl("UK", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Italy", colourVec[index]) == TRUE){
      
      colourVec[index] <- "aquamarine2" 
    } else if(grepl("Spain", colourVec[index]) == TRUE){
      
      colourVec[index] <- "goldenrod3"
    } else if(grepl("France", colourVec[index]) == TRUE){
      
      colourVec[index] <- "royalblue4"
    } else if(grepl("Scotland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("England", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Wales", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Germany", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Netherlands", colourVec[index]) == TRUE){
      
      colourVec[index] <- "orangered"
    } else if(grepl("Czech", colourVec[index]) == TRUE){
      
      colourVec[index] <- "mediumblue"
    } else if(grepl("Greece", colourVec[index]) == TRUE){
      
      colourVec[index] <- "slateblue"
    } else if(grepl("Norway", colourVec[index]) == TRUE){
      
      colourVec[index] <- "purple"
    }else if(grepl("Canada", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deeppink"
    }else if(grepl("Australia", colourVec[index]) == TRUE){
      
      colourVec[index] <- "steelblue3"
    }else if(grepl("USA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    }else if(grepl("Argentina", colourVec[index]) == TRUE){
      
      colourVec[index] <- "lightblue"
    }else if(grepl("India", colourVec[index]) == TRUE){
      
      colourVec[index] <- "lightgreen"
    }else if(grepl("Venezuela", colourVec[index]) == TRUE){
      
      colourVec[index] <- "yellow"
    } else {
      
      colourVec[index] <- "darkgreen"
    }
  }
  
  return(colourVec)
}

