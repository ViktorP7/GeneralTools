#Script to construct phylogenetic tree from CIT data

# Load packages
library(ape)
library(phytools)
library(scales)
library(ggplot2)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaMay2021Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SarahProject/vcfFiles/RAxML_bipartitions.variants"


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

# Root the tree at 505 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 26)

# Drop tips corresponding to low qual and het
todrop <- c(5,10)
newtree <- drop.tip(rootree, todrop)

# Convert branch lengths to SNP values
newtree$edge.length <- round(newtree$edge.length * 1310)

# Find the distances between all isolates
allDist <- cophenetic(newtree)

# Round the distances
allDist <- round(allDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(allDist)){
  
  allDist[index, index] <- NA
}

# Get the herd names
herdNames <- getNames(allDist, "Herd")

# Get current county, birth county names and sameness
countyNames <- getNames(allDist, "CCounty")
birthcountyNames <- getNames(allDist, "BCounty")
sameness <- getNames(allDist, "Same")

# Get VNTR names
vntrNames <- getNames(allDist, "VNTR")

# Make VNTR colours
vntrTips <- makeVNTRCols(vntrNames)

# Simplify the labels
simpleLabels <- deconstructLabels(newtree$tip.label, counties, shortCounties)

# Assign simple labels
newtree$tip.label <- simpleLabels
newtree$tip.label[6] <- "K10"
newtree$tip.label[11] <- "CITP"
newtree$tip.label[12] <- "CIT"

#### Tree plotting (.pdf) ####

# Note different settings to .png

# Save plot as .pdf file (Ireland)
outputFile <- paste("CIT_Tree_27-12-21.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot VNTR tree
plot.phylo(newtree, edge.width = 11, font = 1, label.offset = 0.2, tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 10, no.margin = TRUE)

# Add the SNP scale
add.scale.bar(x=10, y = 5, cex = 12, lwd = 15)
text(x=50, y =5, cex = 12, "SNPs")

dev.off()

#### Growth experiment ####

# Read in tables
pathExp1 <- "C:/Users/UCD/Documents/Lab/CIT Analysis/GrowthExpPass1.csv"
pathExp2 <- "C:/Users/UCD/Documents/Lab/CIT Analysis/GrowthExpPass2.csv"

# Read in table of CVRL isolates
exp1Table <- read.table(pathExp1,
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)
exp2Table <- read.table(pathExp2,
                        header = TRUE,
                        sep = ",",
                        stringsAsFactors=FALSE, 
                        check.names=FALSE)


# Plot the first growth set
colourset=c("red", "blue", "green", "black", "goldenrod3", "violet")
plot(NULL, xlim = c(1,nrow(exp1Table)), ylim = c(0, 1.2), main = "K10 Growth Set 1, 1st Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 2:7){
  lines(exp1Table[,spot], col=colourset[spot-1],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)

plot(NULL, xlim = c(1,nrow(exp1Table)), ylim = c(0, 1.2), main = "CIT Growth Set 1, 1st Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 8:13){
  lines(exp1Table[,spot], col=colourset[spot-7],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)

plot(NULL, xlim = c(1,nrow(exp1Table)), ylim = c(0, 1.2), main = "K10 Growth Set 2, 1st Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 14:19){
  lines(exp1Table[,spot], col=colourset[spot-13],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)

plot(NULL, xlim = c(1,nrow(exp1Table)), ylim = c(0, 1.2), main = "CIT Growth Set 2, 1st Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 20:25){
  lines(exp1Table[,spot], col=colourset[spot-19],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)

plot(NULL, xlim = c(1,5), ylim = c(0, 1.2), main = "K10 Growth, 2nd Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 2:7){
  lines(exp2Table[,spot], col=colourset[spot-1],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)

plot(NULL, xlim = c(1,5), ylim = c(0, 1.2), main = "CIT Growth, 2nd Passage", xlab = "Reading", ylab = "OD600nm")
for(spot in 8:13){
  lines(exp2Table[,spot], col=colourset[spot-7],, type = "b", pch=1)
}
legend("topleft", legend= c("Glucose", "Pyruvate", "Propionate", "Lactate", "Full", "No Carbon"), text.col = colourset)
#### Functions ####

# Function to get the labels names for bryant isolates
getBryantLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  chopVector <- TheTree$tip.label
  
  # Remove trailing _x.vcf.gz
  chopchopVector <- sapply(strsplit(chopVector, split = "_"), function(x) (x[1]))
  
  # Remove the > symbol
  nameVector <- sapply(strsplit(chopchopVector, split = ">"), function(x) (x[2]))
  
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
    if(same == "n/a" || is.na(same) == TRUE || same == "Unknown" || same == "Not available"|| same == "Notavailable"){
      
      yoke <- paste(herd,"*", collapse = NULL)
      
      newtips[index] <- yoke
      
    } else if(same == "None" || same == "Same"){
      
      newtips[index] <- herd 
    }else {
      
      if(birth == "U.K. Import" || birth == "U.K.Import"){
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



