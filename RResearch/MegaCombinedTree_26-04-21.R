# Created 26-04-21

# Load packages
library(ape)
library(phytools)
library(scales)
library(phangorn)
library(seqinr)

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

shortCounties <- c("AnNI", "ArNI", "CW", "CN", "CE", "C", "DL", "DoNI", "D", "FNI", "G", "KY", "KE", "KK", "LS", "LM",
                   "L", "DeNI", "LD", "LH", "MO", "MH", "MN", "OY", "RN", "SO", "T", "TNI", "W", "WH", "WX", "WW")

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaMay2021Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathBigTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/GenieBigMAP/FullRunMay30/RAxML_bipartitions.RaxML-R_02-06-21"
pathNI <- "C:/Users/UCD/Documents/Lab/CVRL MAP/NIMetaOct2020.csv"
pathAus <- "C:/Users/UCD/Documents/Lab/AustraliaMetadata.csv"
pathCan <- "C:/Users/UCD/Documents/Lab/CanadaMetadata.csv"
pathUS <- "C:/Users/UCD/Documents/Lab/USMeta.csv"
pathFasta <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/GenieBigMAP/FullRunMay30/lessallcore.fa"
pathRTT <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/remove-uniformative-sites-master/finalTree-RTT.tsv"

# Read in table of bryant isolates
isoBryantTable <- read.table(pathBryantIso, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of CVRL isolates
isoCVRLTable <- read.table(pathNewIso, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of NI isolates
isoNITable <- read.table(pathNI, header = TRUE, sep = ",", stringsAsFactors=FALSE, check.names=FALSE)

# Read in table of Aus isolates
ausTable <- read.table(pathAus, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

canTable <- read.table(pathCan, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

usTable<- read.table(pathUS, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

# Read in root to tip stuff
rtt <- read.table(pathRTT, header = T, sep = "\t", stringsAsFactors = F, check.names = F)

# Plot rtt stuff
rtt$distance <- rtt$distance*2883
plot(x=rtt$date, y=rtt$distance, main = "Root to Tip Analysis", xlab = "Date", ylab = "SNP Distance")
abline(mod <- lm(rtt$distance ~ rtt$date), col="red")
coef(mod)

# Read in fasta file
treefas <- read.FASTA(pathFasta, type = "DNA")
treenames <- names(treefas)
droppem <- match(BigTree$tip.label[tipstogo], treenames)
droppem=na.omit(droppem)
treefas2 <- treefas[-droppem]
write.FASTA(treefas2,"noSTree.fa")
#### Model Investigation ####

#remove sheep isolates


treefas3 <- as.phyDat(treefas2)

modelResult <- modelTest(treefas3, model = c("JC","HKY","GTR"))

#### Tree file processing ####

# Read in tree
BigTree <- read.tree(pathBigTree)

# Root on another node
BigRoot <- root(BigTree, node = 1301)

# Drop sheep isolates and avium
tipstogo <- c(185:280)
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

# Get us metadata
realNames <- getSRALabels(usTable, noSTree)

noSTree$tip.label <- realNames

# Convert branch lengths to SNP values
noSTree$edge.length <- round(noSTree$edge.length * 31255)

worldDist <- cophenetic(noSTree)
worldDist <- round(worldDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(worldDist)){
  
  worldDist[index, index] <- NA
}

# Make EU colours
worldCols <- newmakeRegionColours(noSTree$tip.label)

# Drop international isolates
dropper <- toDropInternationalTips(noSTree$tip.label)
irishOnlytree <- drop.tip(noSTree, dropper)

simpleLabels <- deconstructLabels(irishOnlytree$tip.label, counties, shortCounties)

# Assign simple labels
irishOnlytree$tip.label <- simpleLabels

irishDist <- cophenetic(irishOnlytree)

# extract releavnt clade
newclade <- extract.clade(irishOnlytree, node = 337)
newcols <- rep("black", 31)
newcols[c(2,3,4,22,31)] <- "red"

plot.phylo(newclade, label.offset = 0.2, no.margin = T, tip.color = newcols)
add.scale.bar(x=80,y=10)


edgelabels(newclade$edge.length)

#extract group f clade
fclade <- extract.clade(irishOnlytree, node=238)
fcols <- rep("black", length(fclade$tip.label))
fcols[30] <- "red"
fcols[c(26,27,28,29,24,20,21,37,38)] <- "blue"

plot.phylo(fclade, label.offset = 0.2, no.margin = T, tip.color = fcols, cex = 0.8)
add.scale.bar(x=125,y=5)
edgelabels(fclade$edge.length)

#extract group e clade
eclade <- extract.clade(irishOnlytree, node=370)
ecols <- rep("black", length(eclade$tip.label))
ecols[c(1,2,12,13)] <- "blue"

plot.phylo(eclade, label.offset = 0.2, no.margin = T, tip.color = ecols, cex = 0.8)
add.scale.bar(x=30,y=7)
edgelabels(eclade$edge.length)
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
outputFile <- paste("World-Tree_28-10-21.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot EU tree
plot.phylo(noSTree, edge.width = 6, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="fan", cex = 30, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","F","G","H"), node = c(1933,1831,1806,1793,1883,1057,1028,1745,1725), frame = "c", cex=10, col = "red")



#Add shaped tip labels
tiplabels(pch = 18, col = worldCols,  cex = 8)

# Add the SNP scale
add.scale.bar(x=-180, y=-250, cex = 10, lwd = 15)
text(x=-150, y=-260, cex = 10, "SNPs")

# Add a legend
legend(x=-250, y=220, legend = c("Ireland", "UK", "Continental Europe","Canada","USA","Australia (& NZ)", "Rest of World"), 
       text.col = c("darkgreen", "firebrick4", "blue", "deeppink","purple",
                    "goldenrod3", "black"), pch = 18,
       col = c("darkgreen", "firebrick4", "blue", "deeppink","purple",
               "goldenrod3", "black"),
       bty = "n", cex = 10, y.intersp = 0.8, title = "Location", title.col = "black")

dev.off()

# Get SNP distances between Ireland and other locations
ie <- grep("darkgreen", worldCols)
aus <- grep("goldenrod3", worldCols)
us <- grep("purple", worldCols)
can <- grep("deeppink", worldCols)
eu <- grep("blue", worldCols)
uk <- grep("firebrick4", worldCols)

ieAus <- getWithinBetween(worldDist[c(ie,aus),c(ie,aus)], worldCols[c(ie,aus)],FALSE) 

ieUS <- getWithinBetween(worldDist[c(ie,us),c(ie,us)], worldCols[c(ie,us)],FALSE) 

ieEU<- getWithinBetween(worldDist[c(ie,eu),c(ie,eu)], worldCols[c(ie,eu)],FALSE) 

ieCan<- getWithinBetween(worldDist[c(ie,can),c(ie,can)], worldCols[c(ie,can)],FALSE) 

ieUK<- getWithinBetween(worldDist[c(ie,uk),c(ie,uk)], worldCols[c(ie,uk)],FALSE) 


boxplot(ieAus$Between, ieUS$Between, ieEU$Between, ieCan$Between, ieUK$Between,
        main = "SNP distances between Irish and non-Irish isolates", 
        names = c("IE/AUS", "IE/US", "IE/EUR", "IE/CAN", "IE/UK"),
        ylab = "SNP Difference")
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
      
      colourVec[index] <- "blue" 
    } else if(grepl("Spain", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("France", colourVec[index]) == TRUE){
      
      colourVec[index] <- "blue"
    } else if(grepl("Scotland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("England", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
    } else if(grepl("Wales", colourVec[index]) == TRUE){
      
      colourVec[index] <- "firebrick4"
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
    }else if(grepl("Canada", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deeppink"
    }else if(grepl("Australia", colourVec[index]) == TRUE){
      
      colourVec[index] <- "goldenrod3"
    }else if(grepl("USA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "purple"
    }else if(grepl("Argentina", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    }else if(grepl("India", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    }else if(grepl("Venezuela", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    }else if(grepl("New Zealand", colourVec[index]) == TRUE){
      
      colourVec[index] <- "goldenrod3"
    } else {
      
      colourVec[index] <- "darkgreen"
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
    } else if(grepl("SRR", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    } else if(grepl("Australia", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    }
  }
  return(dropVector)
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

# Function to pull out within and between SNP distances
getWithinBetween <- function(mat, name, shuffle){
  
  if(shuffle == TRUE){
    
    shuffler <- sample(name, replace = FALSE) 
  } else {
    
    shuffler <- name
  }
  
  # Create list to store info in
  withinBetween <- list()
  
  # Create keys
  keys <- c("CellLabel", "Within", "Between")
  
  # Loop thru rows of input matrix
  for(row in 1:nrow(mat)){
    
    # Loop thru columns
    for(col in 1:ncol(mat)){
      
      # Check if value present in cell
      if(is.na(mat[row,col]) == TRUE){
        
        next
      } else {
        
        # Are they the same VNTR
        if(shuffler[row] == shuffler[col]){
          
          # Get the cell and status
          status <- paste(row, col, "Within", sep = "_")
          
          # Store cell location and status
          withinBetween[[keys[1]]] <- append(withinBetween[[keys[1]]], status)
          
          # Store distance value
          withinBetween[[keys[2]]] <- append(withinBetween[[keys[2]]], mat[row,col])
          
        }else{
          
          # Get the cell and status
          status <- paste(row, col, "Between", sep = "_")
          
          # Store cell location and status
          withinBetween[[keys[1]]] <- append(withinBetween[[keys[1]]], status)
          
          # Store distance value
          withinBetween[[keys[3]]] <- append(withinBetween[[keys[3]]], mat[row,col])
        }
        
      }
      
    }
    
  }
  return(withinBetween)        
}