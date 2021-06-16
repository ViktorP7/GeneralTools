# Created 01-06-21

# Load packages
library(ape)
library(phytools)
library(scales)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaMay2021Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathC10Tree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Cork10Fastqs/C10genie/RAxML_bipartitions.RaxML-R_01-06-21"

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
c10Tree <- read.tree(pathC10Tree)

# Re-root
c10reRoot <- root(c10Tree, node = 105)

# Get the Bryant names
c10realNames <- altinputgetBryantLabels(isoBryantTable, c10reRoot)

# Update the names in the tree
c10reRoot$tip.label <- c10realNames

# Get metadata for CVRL isolates
c10realNames <- altgetCVRLLabels(isoCVRLTable, c10reRoot)

# Update the names in the tree
c10reRoot$tip.label <- c10realNames

#Get rid of non-relevant tips
#dropem <- c(17,19,20,95,96,135,136,183,184)
#onlytree <- drop.tip(irishOnlytree, dropem)

# Convert branch lengths to SNP values
#onlytree$edge.length <- round(onlytree$edge.length * 7843)
c10reRoot$edge.length <- round(c10reRoot$edge.length * 135)

# Find the distances between all isolates
#allDist <- cophenetic(onlytree)
c10Dist <- cophenetic(c10reRoot)

# Round the distances
c10Dist <- round(c10Dist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(c10Dist)){
  
  c10Dist[index, index] <- NA
}

# Get the herd names
c10herdNames <- getNames(c10Dist, "Herd")

# Get current county, birth county names and sameness
c10countyNames <- getNames(c10Dist, "CCounty")
c10birthcountyNames <- getNames(c10Dist, "BCounty")
c10sameness <- getNames(c10Dist, "Same")

# Get VNTR names
c10vntrNames <- getNames(c10Dist, "VNTR")

# Make VNTR colours
c10vntrTips <- makeVNTRCols(c10vntrNames)

# Simplify the labels
c10simpleMat <- altdeconstructLabels(c10reRoot$tip.label, counties, shortCounties)

# Assign simple labels
c10reRoot$tip.label <- c10simpleMat[1,]

#### Tree plotting (.png) ####

# Save plot as .png file (Ireland)
outputFile <- paste("VNTR_Tree_10-05-21.png", sep="")
png(outputFile, height=10000, width=6000)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 11, font = 1, label.offset = 0.2, 
           tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 5, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(317,285,267,259,249,222,213,199), frame = "n", cex=15, adj = c(1,1), col = "red")


# Add the SNP scale
add.scale.bar(x=10, y = 5, cex = 8, lwd = 15)
text(x=40, y =5, cex = 8, "SNPs")

# Add a legend
legend(x=9, y=160, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types", title.col = "black")

dev.off()

#### Tree plotting (.pdf) ####

# Save plot as .pdf file (Ireland)
outputFile <- paste("C10_Tree_02-06-21.pdf", sep="")
pdf(outputFile, height=50, width=50)

# Plot VNTR tree
plot.phylo(c10reRoot, edge.width = 17, font = 2, label.offset = 0.5, tip.color = c10vntrTips,
           align.tip.label = TRUE, type="phylogram", cex = 3.5, no.margin = TRUE)

tiplabels(pch = 18, frame = "n", col = c10simpleMat[2,], cex=10)
# Add the SNP scale
add.scale.bar(x=2, y = 25, cex = 10, lwd = 17)
text(x=10,y=25, cex = 10, "SNPs")
text(x=5,y=15, cex=5, "*-* Unique Animal")
text(x=5,y=12, cex=5, "*A* Serial Animal")


# Add a legend
legend(x=0, y=100, legend = c("(42332228) - 1", "(32332228) - 2"), 
       text.col = c("red", "deepskyblue3"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types", title.col = "black")
legend(x=2, y=75, legend = c("2014", "2016", "2017", "2018", "2019"), 
       text.col = c( "goldenrod3", "steelblue3", "palegreen3","gray40", "orchid"),
       pch = c(23,23,23,23,23), pt.bg = c("goldenrod3",  "steelblue3", "palegreen3","gray40", "orchid"), pt.cex = 15,
       bty = "n", cex = 10, y.intersp = 0.8, title = "Isolation Year", title.col = "black")

dev.off()

#### Functions ####

# Function to get the labels names for bryant isolates
altinputgetBryantLabels <- function(isoTable, TheTree){
  
  # Create a vector to match the tip labels vector in the tree
  nameVector <- TheTree$tip.label
  
  # Loop thru the table
  for(row in 1:nrow(isoTable)){
    
    # Loop thru the name vector
    for(index in 1:length(nameVector)){
      
      # Check if the current accession is present in the big table
      if(isoTable[row,"Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row,"Country of origin"], "_", splitter)
        nameVector[index] <- newname
      }else if(isoTable[row,"Secondary Accession"] == nameVector[index]){
        
        splitter <- strsplit(isoTable[row,"Reference No"], split = "I")[[1]][2]
        newname <- paste(isoTable[row,"Country of origin"], "_", splitter)
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
altgetCVRLLabels <- function(isoTable, TheTree){
  
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
                         isoTable[row,"INMV Group"], "_", isoTable[row, "Birth County"], "_", isoTable[row, "Number of moves to herd of sampling"],"_",isoTable[row,"AnimalNum"])
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

# Function to simplify the labels
altdeconstructLabels <- function(tiplabel, counties, shortCounties){
  
  # Create output frame
  outmat <- matrix(nrow = 2, ncol = length(tiplabel))
  
  
  # Loop thru the tips and cut them down
  for(index in 1:length(tiplabel)){
    
    if(grepl("Scotland", tiplabel[index]) == TRUE || grepl("ref", tiplabel[index]) == TRUE){
      
      outmat[1,index] <- tiplabel[index]
      outmat[2,index] <- "darkred"
    } else{
    
      # Split up the different parts of the tip label
      one <- strsplit(tiplabel[index], split = "_")[[1]][2]
      two <- strsplit(tiplabel[index], split = "_")[[1]][3]
      notdate <- strsplit(tiplabel[index], split = "_")[[1]][1]
      date <- strsplit(notdate, split = "-")[[1]][1]
      animal <- strsplit(tiplabel[index],split="_")[[1]][7]
    
      # Find which index in counties the tip county is and store shortened version
      short <- shortCounties[which(one == counties)]
    
      # Store herd
      herd <- paste(short,two)
    
      # Store birth location
      birth <- strsplit(tiplabel[index], split = "_")[[1]][5]
    
    
      # Store sameness
      same <- strsplit(tiplabel[index], split = "_")[[1]][6]
    
      # Check if it's the same county
      if(same == "n/a" || is.na(same) == TRUE || same == "Unknown" || same == "Not available"|| same == "Notavailable"){
      
        yoke <- paste(notdate,herd,"*",animal, collapse = NULL)
      
        outmat[1,index] <- yoke
      
      } else if(same == "None" || same == "Same"){
      
        outmat[1,index] <- paste(notdate, herd, animal) 
      }else {
      
        if(birth == "U.K. Import" || birth == "U.K.Import"){
          shortB <- "UK"
        } else{
        
          shortB <- shortCounties[which(birth == counties)]
        }
        birthstring <- paste("(",shortB,")", sep = "")
      
        thing <- paste(notdate, herd, animal, birthstring )
      
        outmat[1,index] <- thing
      }
    
      if(date == "13"){
      
        co = "firebrick3"
      } else if(date == "14"){
      
        co = "goldenrod3"
      } else if(date == "18"){
      
        co = "gray40"
      } else if(date == "16"){
      
        co = "steelblue3"
      } else if(date == "17"){
      
        co = "palegreen3"
      } else if (date == "19"){
      
        co = "orchid"
      }
    
      outmat[2,index] <- co
    }  
  }
  
  return(outmat)
}

