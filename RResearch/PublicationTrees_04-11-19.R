# 04-11-19 Created script

#### Package import and file loading ####

# Load packages
library(ape)
library(phytools)
library(scales)
library(sp)
library(raster)
library(gplots)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MAP-Metadata-Formatted-May19.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_30-09-19"


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

# Root the tree at 473 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 491)

# Drop tips for far away ancestors 
dropNumbers <- c(465,466)
droppedTree <- drop.tip(rootree, dropNumbers)

# Extract the clade that doesn't have all the distant sheep and cows
extractedTree <- extract.clade(droppedTree, node = 507)

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(extractedTree$tip.label)
irishOnlytree <- drop.tip(extractedTree, dropper)

# Get rid of non-relevant tips
dropem <- c(174,175,132,133,93,94,95,96,19,2)
onlytree <- drop.tip(irishOnlytree, dropem)

# Correct mislabel
onlytree$tip.label[53] <- "17-5652_Cork_10_1_Cork"

# Convert branch lengths to SNP values
onlytree$edge.length <- onlytree$edge.length * 68928

# Get the rounded values o the lengths
roundedSNPs <- round(onlytree$edge.length)

# Assign floored SNPs
onlytree$edge.length <- roundedSNPs

# Remove K10
useTree <- drop.tip(onlytree, 2)

# Find the distances between all isolates
allDist <- cophenetic(useTree)

# Round the distances
allDist <- round(allDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(allDist)){
  
  allDist[index, index] <- NA
}

# Get the herd names
herdNames <- getNames(allDist, "Herd")

# Get rid of all Cork 10 isolates except two
notTheseTips <- cork10Subsampler(herdNames, useTree$tip.label)

# Drop'em
subbedTree <- drop.tip(useTree, notTheseTips)

# Find the distances between all isolates
newDist <- cophenetic(subbedTree)

# Round the distances
newDist <- round(newDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(newDist)){
  
  newDist[index, index] <- NA
}

# Get VNTR names
vntrNames <- getNames(newDist, "VNTR")

# Get herd names
newHerdNames <- getNames(newDist, "Herd")

# Make VNTR colours
vntrTips <- makeVNTRCols(vntrNames)

# Simplify the labels
simpleLabels <- deconstructLabels(subbedTree$tip.label)

# Assign simple labels
subbedTree$tip.label <- simpleLabels

#### VNTR tree plotting ####

# Save plot as .pdf file
outputFile <- paste("VNTR_Tree.pdf", sep="")
pdf(outputFile, height=35, width=20)

# Set margins to nothing
currentMar <- par()$mar
par(mar=c(0,0,0,0))

# Plot VNTR tree
plot.phylo(subbedTree, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 1.6)

# Add the SNP scale
add.scale.bar(x=100, y = 2, cex = 4)
text(x=135, y =2, cex = 4, "SNPs")

# Add a legend
legend(x=115, y=110, legend = c("1", "2", "3", "13", "116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 4, y.intersp = 1)
text(x=130, y=110, "INMV Type", cex = 4.5)


# Reset the margins
par(mar=currentMar)

dev.off()

#### Processing of files for map building ####
# Set path variables
pathCoords <- "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Get the numbers of isolates
sampleNumbers <- numberSamples(newHerdNames)

# Get polygon coordinates for map plotting
polygonCoords <- getSpatialPolygons(counties, pathCoords)

# Calculate limits of plot
ranges <- mapSpatialLimits(polygonCoords, counties)

# Calculate the max number of samples
maxNSamples <- sampleMax(sampleNumbers)

# Save plot as .pdf file
outputFile <- paste("IrishMap.pdf", sep="")
pdf(outputFile, height=20, width=15)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Add legend
legend("topleft", legend = c("1","2","3","4","5","6"), title="Isolates per Herd", bty="n", cex=3,
       pch = c(21,21,21,21,21,21), col = alpha("red", 0.9), pt.bg = alpha("red", 0.65), 
       pt.cex = c(2,3,4,5,6,7))

# Plot county polygons and related sample data
polygonsSpatialData(polygonCoords, sampleNumbers, counties, newHerdNames)

dev.off()

#### Temporal Plot ####

# Get the count matrix
counterMatrix <- getYearCounts(newHerdNames, newDist)

# Plot the figure
heatmap.2(counterMatrix, dendrogram = "none", trace = "none", Rowv = FALSE,
          Colv = FALSE, col = )

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
      }else if(nameVector[index] == "14-2662"){
        
        nameVector[index] <- "14-2622"
      }else if(nameVector[index] == "17-5652"){
        
        nameVector[index] <- "17-6652"
        
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
        
        herd <- strsplit(isoTable[row,"Herd Ref"], split = " ")[[1]][2]
        
        newname <- paste(nameVector[index], "_", isoTable[row, "Herd Location"], "_", herd, "_",
                         isoTable[row,"INMV Group"], "_", isoTable[row, "County of Birth"])
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

# Function to remove Cork 10 apart from 2 random isolates
cork10Subsampler <- function(herdNames, tips){
  
  # Vector to store indexes
  indexes <- c()
  
  # Loop thru each herd
  for(herd in unique(herdNames)){
    
    # Check if the herd is Cork 10
    if(herd == "Cork _ 10"){
      
      # Subsample from the big herd
      subs <- sample(which(herdNames %in% herd), 2, replace = FALSE)
      
      # Place subsample indexes in
      indexes <- append(indexes, subs)
    } else{
      
      # Pull out the indexes and place in indexes vector
      indexes <- append(indexes, which(herdNames %in% herd))
    }
  }
  
  # Create vector of indexes from tips 
  tipper <- c(1:length(tips))
  
  # Get rid of indexes from tips
  tipsOut <- tipper[-indexes]
  
  return(tipsOut)
  
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
  }
}

# Function to create tip labels colours based on VNTR
makeVNTRCols <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    if(is.na(colourVec[index]) == TRUE){
      
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
deconstructLabels <- function(tiplabel){
  
  # Copy vector
  newtips <- rep(NA, length(tiplabel))
  
  # Loop thru the tips and cut them down
  for(index in 1:length(tiplabel)){
    
    # Split up the different parts of the tip label
    one <- strsplit(tiplabel[index], split = "_")[[1]][2]
    two <- strsplit(tiplabel[index], split = "_")[[1]][3]
    
    # Store herd
    herd <- paste(one,two)
    
    # Store birth location
    birth <- strsplit(tiplabel[index], split = "_")[[1]][5]
    
    # Check if the birth county is equal to current
    if(grepl(one, birth) == TRUE){
      
      newtips[index] <- herd
      
    }else if(grepl("n/a", birth) == TRUE || is.na(birth) == TRUE){
      
      yoke <- paste(herd,"*", collapse = NULL)
      
      newtips[index] <- yoke
    } else {
      
      thing <- paste(herd,"(",birth,")")
      
      newtips[index] <- thing
    }
    
  }
  
  return(newtips)
}

# Function to count amount of isolates and herds per each county
numberSamples <- function(isolates){ 
  
  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()
  countyHerds <- list()
  
  # Split off the county name away from herd number
  for(index in 1:length(isolates)){ 
    
    # Split away the herd number
    countyName <- strsplit(isolates[index], split = " ")[[1]][1]
    herdSampled <- strsplit(isolates[index], split = " ")[[1]][3]
    
    # Check if we have encountered the current county before
    if(is.null(countySamples[[countyName]]) == TRUE){
      
      # add it into the list and count once
      countySamples[[countyName]] <- c(1, 1) 
      countyHerds[[countyName]] <- c(herdSampled)
      
      # if already there add an extra 1 to the count for each time seen
    } else { 
      
      countySamples[[countyName]][1] <- countySamples[[countyName]][1] + 1
      
      if(herdSampled %in% countyHerds[[countyName]] == FALSE){
        countyHerds[[countyName]] <- c(countyHerds[[countyName]], herdSampled)
        countySamples[[countyName]][2] <- length(countyHerds[[countyName]])
      }
    }
  }
  return(countySamples)
}

# Function to get spatial polygons
getSpatialPolygons <- function(counties, path) { # input is county list
  
  polygonCoords <- list() # empty list to store gps data for each county
  
  
  for(index in 1:length(counties)) { # for each county
    
    fileName <- paste(path, "PolygonCoords_", counties[index], ".txt", sep="" )
    # paste general path to counties, and the .txt ending to get path for each
    
    # read in the values for each county into the appropriate list segment
    polygonCoords[[counties[index]]] <- SpatialPolygons(list(Polygons(list(Polygon(read.table(fileName, 
                                                   header = TRUE, sep = "\t"))), "x")))
  }
  return(polygonCoords)
}

# Function to calculate map ranges
mapSpatialLimits <- function(polygonCoords, counties) {
  
  # Initialise vectors to store the mins and maxes of each county
  minX <- rep(NA, length(counties))
  maxX <- rep(NA, length(counties))
  minY <- rep(NA, length(counties))
  maxY <- rep(NA, length(counties))
  
  # For each county, calculate the mins and maxes, and populate the empty vectors
  for(index in 1:length(counties)) {
    
    minX[index] <- extent(polygonCoords[[counties[index]]])[1]
    maxX[index] <- extent(polygonCoords[[counties[index]]])[2]
    minY[index] <- extent(polygonCoords[[counties[index]]])[3]
    maxY[index] <- extent(polygonCoords[[counties[index]]])[4]
  }
  
  # Get the overall mins and maxes of X and Y to be able to set plot limits
  ranges <- c(min(minX), max(maxX), min(minY), max(maxY))
  
  return(ranges)
}

# Calculate max samples
sampleMax <- function(sampleNumbers) {
  currentmaxValue <- 1 # set current max = 1
  
  keys <- names(sampleNumbers) # assign names as keys to access list
  
  for(index in 1:length(keys)) { # for every element in the list
    if(sampleNumbers[[keys[index]]][1] > currentmaxValue) { # if they are greater than 1
      currentmaxValue = sampleNumbers[[keys[index]]][1] # current max is changed to that
    }
  }
  
  return(currentmaxValue)
}

# Fill polygons on map plot
polygonsSpatialData <- function(polygonCoords, sampleNumbers, counties, herdnames) {
  
  for(index in 1:length(counties)) {
 
    # Add county name
    if(is.null(sampleNumbers[[counties[index]]]) == TRUE ){
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
              border = "black", add = TRUE)

    }else{ # counties that are present on the sample list
      
      # Get number of samples for current county
      nSamples <- sampleNumbers[[counties[index]]][1]
      nHerdsSamples <- sampleNumbers[[counties[index]]][2]
      
      # Plot the polygon of the current county
      plot(polygonCoords[[counties[index]]],
              border = "black", add = TRUE)
      
      # Create a vector of herds in the current county
      currentHerds <- unique(herdnames[grep(counties[index], newHerdNames)])
      
      # For each herd in the current county herds, count how many isolates and plot point
      for(herd in currentHerds){
        
        pointSize <- length(grep(herd, herdnames))
        
        points(spsample(polygonCoords[[counties[index]]], n=1, type = "random"),
               pch = 21, col = alpha("red", 0.9), bg = alpha("red", 0.65), cex = 1+pointSize)
      }
    }
  }
}

# Function to get a table of counts of isolates for each year
getYearCounts <- function(herdnames, matrix){
  
  # Create a new vector of length herdnames
  yearVector <- rep(NA, length(herdnames))
  
  # Loop thru columns and extract years
  for(col in 1:ncol(matrix)){
    
    yearVector[col] <- strsplit(colnames(matrix)[col], split = "-")[[1]][1]
  }
  
  # Paste the date into the herd vector
  comboVector <- paste(yearVector,"-",herdnames)
  
  # Tabulate the counts of each
  herdYears <- table(comboVector)
  
  # Make a matrix of length herds and width years
  countMatrix <- matrix(data = 0, nrow = length(unique(herdnames)), ncol = 5, 
                        dimnames = list(sort(unique(herdnames)), c("2012","2013","2014","2016","2017")))
  
  # Get the name of the table
  herdos <- names(herdYears)
  
  # Loop thru the names
  for(index in 1:length(herdos)){
    
    # Split off the year and herd
    year <- strsplit(herdos[index], split = "-")[[1]][1]
    place <- strsplit(herdos[index], split = "-")[[1]][2]
    
    # Match column
    currentCol <- grep(trimws(year), colnames(countMatrix))
    
    # Match row
    currentRow <- grep(trimws(place), rownames(countMatrix))
    
    # Add a count to the relevant matrix slot
    countMatrix[currentRow, currentCol] <- herdYears[index]
    
  }
  
  return(countMatrix)
}