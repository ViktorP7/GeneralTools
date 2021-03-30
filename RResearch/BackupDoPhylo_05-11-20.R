	
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
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_22-10-20"
pathNI <- "C:/Users/UCD/Documents/Lab/CVRL MAP/NIMetaOct2020.csv"

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
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)

# Read in table of NI isolates
isoNITable <- read.table(pathNI,
                         header = TRUE,
                         sep = ",",
                         stringsAsFactors=FALSE, 
                         check.names=FALSE)

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

shortCounties <- c("AnNI", "ArNI", "CW", "CN", "CE", "C", "DL", "DoNI", "D", "FNI", "G", "KY", "KE", "KK", "LS", "LM",
                   "L", "DeNI", "LD", "LH", "MO", "MH", "MN", "OY", "RN", "SO", "T", "TNI", "W", "WH", "WX", "WW")

#### Tree file processing ####

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

# Get metadata for NI isolates
realNames <- getNILabels(isoNITable, TheTree)

# Update names
TheTree$tip.label <- realNames

# Root the tree at 505 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 750)

# Drop tips for far away ancestors 
dropNumbers <- c(351:356)
droppedTree <- drop.tip(rootree, dropNumbers)

# Extract the clade that doesn't have all the distant sheep and cows
extractedTree <- extract.clade(droppedTree, node = 478)

# Get rid of non-EU isolates for international tree
dropInternational <- toDropNonEUTips(extractedTree$tip.label)
euOnlyTree <- drop.tip(extractedTree, dropInternational)

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(extractedTree$tip.label)
irishOnlytree <- drop.tip(extractedTree, dropper)

# Get rid of non-relevant tips
dropem <- c(1,3,111,113,135,182,183,200,201)
onlytree <- drop.tip(irishOnlytree, dropem)

# Convert branch lengths to SNP values
onlytree$edge.length <- round(onlytree$edge.length * 57086)
euOnlyTree$edge.length <- round(euOnlyTree$edge.length * 57086)

# Find the distances between all isolates
allDist <- cophenetic(onlytree)
euDist <- cophenetic(euOnlyTree)

# Round the distances
allDist <- round(allDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(allDist)){
  
  allDist[index, index] <- NA
}

# Get the herd names
herdNames <- getNames(allDist, "Herd")
euHerds <- getNames(euDist, "Herd")

# Get current county, birth county names and sameness
countyNames <- getNames(allDist, "CCounty")
birthcountyNames <- getNames(allDist, "BCounty")
sameness <- getNames(allDist, "Same")

# Get VNTR names
vntrNames <- getNames(allDist, "VNTR")

# Make VNTR colours
vntrTips <- makeVNTRCols(vntrNames)

# Make EU colours
euCols <- makeRegionColours(euOnlyTree$tip.label)

# Simplify the labels
simpleLabels <- deconstructLabels(onlytree$tip.label, counties, shortCounties)

# Deal with stubborn label
simpleLabels[227] <- "14-3960 DL 1 (UK)"
simpleLabels[149] <- "14-3110 TNI CaA *"
simpleLabels[220] <- "19-8936 G 2 (CE)"

# Assign simple labels
onlytree$tip.label <- simpleLabels

#### Tree plotting (.png) ####

# Save plot as .png file (Ireland)
outputFile <- paste("VNTR_Tree_20-10-20.png", sep="")
png(outputFile, height=5000, width=4500)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 15, font = 1, label.offset = 0.2, 
           tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 3.8, no.margin = TRUE)

# Add the SNP scale
add.scale.bar(x=0, y = -1, cex = 8, lwd = 17)
text(x=30, y =-1, cex = 8, "SNPs")

# Add a legend
legend(x=7, y=122, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 8, y.intersp = 0.7, title = "INMV Types")

#segments(x0=130, y0=c(124,118,73.5,68,54,36.5,29,9.5), y1=c(118.5,74.5,69,55,38,30,11,2), lwd = 5)
#text(x=132, y=c(120.5,95,70.5,61,45,33,20,6), c("A","B","C","D","E","F","G","H"), cex=8)

dev.off()

# Make EU plot
outputFile <- paste("EU-Tree_20-10-20.png", sep="")
png(outputFile, height=5000, width=4500)

# Plot EU tree
plot.phylo(euOnlyTree, edge.width = 10, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)

#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 10)

# Add the SNP scale
add.scale.bar(x=0, y=0, cex = 8, lwd = 15)
text(x=31, y=0, cex = 8, "SNPs")

# Add a legend
legend(x=120, y=220, legend = c("Ireland", "UK", "England", "Scotland", "Wales",
                                "Italy", "Spain", "France", "Germany", "Netherlands",
                                "Czech Rep.", "Greece", "Norway"), 
       text.col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
                    "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
                    "mediumblue", "slateblue", "purple"), pch = 18,
       col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
               "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
               "mediumblue", "slateblue", "purple"),
       bty = "n", cex = 8.8, y.intersp = 0.8, title = "Country")

# Add the group idents
#segments(x0=154, y0=c(213.5,202.5,142,111.5,93.5,57.5,48,18), y1=c(204.5,145,114,97.5,59,49.5,19.5,6.5), lwd = 5)
#text(x=156, y=c(209,170,129,104,75,53,34,12), c("A","B","C","D","E","F","G","H"), cex=8)

dev.off()


#### Tree plotting (.pdf) ####

# Note different settings to .png

# Save plot as .pdf file (Ireland)
outputFile <- paste("VNTR_Tree_22-10-20.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 11, font = 1, label.offset = 0.2, tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 2, no.margin = TRUE)

# Add the SNP scale
add.scale.bar(x=10, y = 10, cex = 8, lwd = 15)
text(x=40, y =10, cex = 8, "SNPs")

# Add a legend
legend(x=9, y=220, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.7, title = "INMV Types", title.col = "black")

# Add the group idents
#segments(x0=120, y0=c(124,118,73.5,68,54,36.5,29,9.5), y1=c(118.5,74.5,69,55,38,30,11,2), lwd = 5)
#text(x=122, y=c(120.5,95,70.5,61,45,33,20,6), c("A","B","C","D","E","F","G","H"), cex=8)

dev.off()

# Make EU plot pdf
outputFile <- paste("EU-Tree_22-10-20.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot EU tree
plot.phylo(euOnlyTree, edge.width = 11, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)

#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 9)

# Add the SNP scale
add.scale.bar(x=0, y=0, cex = 8, lwd = 15)
text(x=31, y=0, cex = 8, "SNPs")

# Add a legend
legend(x=120, y=320, legend = c("Ireland", "UK", "England", "Scotland", "Wales",
                                "Italy", "Spain", "France", "Germany", "Netherlands",
                                "Czech Rep.", "Greece", "Norway"), 
       text.col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
                    "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
                    "mediumblue", "slateblue", "purple"), pch = 18,
       col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
               "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
               "mediumblue", "slateblue", "purple"),
       bty = "n", cex = 8.8, y.intersp = 0.8, title = "Country", title.col = "black")

# Add the group idents
#segments(x0=153.5, y0=c(213.5,202.5,142,111.5,93.5,57.5,48,18), y1=c(204.5,145,114,97.5,59,49.5,19.5,6.5), lwd = 5)
#text(x=154.5, y=c(209,170,129,104,75,53,34,12), c("A","B","C","D","E","F","G","H"), cex=8)

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
    } else if(grepl("SRR", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    }
  }
  return(dropVector)
}

# Function to create vector with non-EU tips to drop
toDropNonEUTips <- function(tiplabel){
  
  # Create vector to store index values of tips to be dropped
  dropVector <- c()
  
  # Loop thru tip labels and drop as required
  for(index in 1:length(tiplabel)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("NewZealand", tiplabel[index]) == TRUE){
      
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
    } else if(grepl("SRR", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    }
  }
  return(dropVector)
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

# Function to create tip labels colours based on VNTR
makeVNTRCols <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    if(is.na(colourVec[index]) == TRUE || colourVec[index] == "n/a"){
      
      colourVec[index] <- "grey30"
    } else if(colourVec[index] == "1"){
      
      colourVec[index] <- "red"
    } else if(colourVec[index] == "2"){
      colourVec[index] <- "deepskyblue3"
    } else if(colourVec[index] == "3"){
      colourVec[index] <- "darkorange3"
    } else if(colourVec[index] == "13"){
      colourVec[index] <- "black"
    } else if(colourVec[index] == "116"){
      colourVec[index] <- "darkgreen"
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
      
      colourVec[index] <- "steelblue3"
    } else if(grepl("England", colourVec[index]) == TRUE){
      
      colourVec[index] <- "lightpink2"
    } else if(grepl("Wales", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deeppink"
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
    } else {
      
      colourVec[index] <- "darkgreen"
    }
  }
  
  return(colourVec)
}

# Function to simplify the labels
deconstructLabels <- function(tiplabel, counties, shortCounties){
  
  # Copy vector
  newtips <- rep(NA, length(tiplabel))
  
  # Loop thru the tips and cut them down
  for(index in 1:length(tiplabel)){
    
    # Split up the different parts of the tip label
    one <- strsplit(tiplabel[index], split = "_")[[1]][2]
    two <- strsplit(tiplabel[index], split = "_")[[1]][3]
    date <- strsplit(tiplabel[index], split = "_")[[1]][1]
    
    # Find which index in counties the tip county is and store shortened version
    short <- shortCounties[which(one == counties)]
    
    # Store herd
    herd <- paste(date,short,two)
    
    # Store birth location
    birth <- strsplit(tiplabel[index], split = "_")[[1]][5]
    
    
    # Store sameness
    same <- strsplit(tiplabel[index], split = "_")[[1]][6]
    
    # Check if it's the same county
    if(same == "n/a" || is.na(same) == TRUE || same == "Unknown" || same == "Not available"){
      
      yoke <- paste(herd,"*", collapse = NULL)
      
      newtips[index] <- yoke
      
    } else if(same == "None"){
      
      newtips[index] <- herd 
    }else {
      
      if(birth == "U.K. Import"){
        shortB <- "UK"
      } else{
        
        shortB <- shortCounties[which(birth == counties)]
      }
      birthstring <- paste("(",shortB,")", sep = "")
      
      thing <- paste(herd, birthstring)
      
      newtips[index] <- thing
    }
    
  }
  
  return(newtips)
}

