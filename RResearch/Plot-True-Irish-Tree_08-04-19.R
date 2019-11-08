#### Changelog ####

# 11-03-19 modified for a RAxML tree input file
# 25-03-19 modified for extra isolates
# 08-04-19 modified for better root and to find VNTR group distances
# 07-05-19 modified for new metadata
# 12-06-19 modified for new data
# 26-06-19 cleaned up code
# 24-07-19 cleaned up code & updated with mapping functions
# 07-08-19 modified for new data and reworked getWB function to make it faster
# 01-10-19 modified for new data, added VNTR plotting function

#### Loading of packages & files ####

# Load packages
library(ape)
library(phytools)
library(scales)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MAP-Metadata-Formatted-May19.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_30-09-19"
pathCoords <- "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")


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
#### Tree-vis ####

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

# Plot the tree for visualisation to determine clade extraction
plot.phylo(TheTree, edge.width = 0.2, font = 1, 
           show.tip.label = TRUE, cex = 0.1)
nodelabels(cex = 0.05, frame = "none")

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

# Convert branch lengths to SNP values
onlytree$edge.length <- onlytree$edge.length * 68928

# Get the rounded values o the lengths
roundedSNPs <- round(onlytree$edge.length)

# Assign floored SNPs
onlytree$edge.length <- roundedSNPs

# Correct mislabel
onlytree$tip.label[53] <- "17-5652_Cork_10_1_Cork"

# Get the colours
tipColours <- makeIrishRegionColours(onlytree$tip.label)

# Plot Irish tree
plotIrishTree(onlytree, tipColours)

# Get polygon coordinates for map plotting
polygonCoords <- getPolygonCoords(counties, pathCoords)

# Calculate limits of plot
ranges <- mapLimits(polygonCoords, counties)

# Plot Irish fan tree
plotIrishFan(onlytree, tipColours, polygonCoords, counties, ranges)

# Plot unrooted tree
plotIrishUnrooted(onlytree, tipColours, polygonCoords, counties, ranges)

# Extract Cork 10 herd only
cork10 <- extract.clade(onlytree, node = 190)

# Plot Cork problem herd
plotProblemHerd(cork10)


#### Post Tree Analysis ####

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

# Get names of matrix
nameVec <- getNames(allDist, "VNTR")

# Get the within and between distances
distList <- getWithinBetween(allDist, nameVec, FALSE)

# Plot the within between plot
plotWB(distList, "VNTR types")

# Run permutation and plot the plot
runplotPer(distList, "VNTR", allDist, nameVec, TRUE, 10000, 5)

#### Do with herd names now instead of VNTR####
# Get the herd names
herdNames <- getNames(allDist, "Herd")

# Get the within and between distances for herds
distHerdList <- getWithinBetween(allDist, herdNames, FALSE)

# Plot the within between plot
plotWB(distHerdList, "herds")

# Run permutation and plot the plot
runplotPer(distHerdList, "Herd", allDist, herdNames, TRUE, 10000, 30)

#### Now with county ####
countyNames <- getNames(allDist, "CCounty")

# Get the within and between distances for counties
distCountyList <- getWithinBetween(allDist, countyNames, FALSE)

# Plot the within between plot
plotWB(distCountyList, "counties")

# Run permutation and plot the plot
runplotPer(distCountyList, "County", allDist, countyNames, TRUE, 10000, 30)

#### Look at the birth county now ####
birthcountyNames <- getNames(allDist, "BCounty")

# Change a few manually
birthcountyNames <- birthcountyNames[-c(67, 104, 142)]

# Edit allDist to exclude the excluded isolates
newallDist <- allDist[-c(67,104,142),-c(67,104,142)]

# Get the within and between distances for counties
distbirthCountyList <- getWithinBetween(newallDist, birthcountyNames, FALSE)

# Plot the within between plot
plotWB(distbirthCountyList, "birth counties")

# Run permutation and plot the plot
runplotPer(distbirthCountyList, "Birth County", newallDist, birthcountyNames, TRUE, 10000, 30)

#### Make a VNTR tree ####
vntrNames <- getNames(allDist, "VNTR")

vntrTips <- makeVNTRCols(vntrNames)

plotVNTRTrees(onlytree, vntrTips)

#### Perform subsampling of Ireland ####

# Find tips to drop
notTheseTips <- herdSubsampler(herdNames, useTree$tip.label)

# Drop'em
subbedTree <- drop.tip(useTree, notTheseTips)

# Get rid of the colours too
newCols <- makeIrishRegionColours(subbedTree$tip.label)

# Plot the subbed tree
plotIrishTree(subbedTree, newCols)
plotIrishFan(subbedTree, newCols, polygonCoords, counties, ranges)

# Check genetic distances
subDist <- cophenetic(subbedTree)

# Round the distances
subDist <- round(subDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(subDist)){
  
  subDist[index, index] <- NA
}

# Get names of matrix
subVNTR <- getNames(subDist, "VNTR")

# Get the within and between distances
subVNTRList <- getWithinBetween(subDist, subVNTR, FALSE)

# Plot the within between plot
plotWB(subVNTRList, "VNTR types")

subHerds <- getNames(subDist, "Herd")

# Get the within and between distances for herds
subHerdList <- getWithinBetween(subDist, subHerds, FALSE)

# Plot the within between plot
plotWB(subHerdList, "herds")

subCounties <- getNames(subDist, "CCounty")

# Get the within and between distances for counties
subCountyList <- getWithinBetween(subDist, subCounties, FALSE)

# Plot the within between plot
plotWB(subCountyList, "counties")

subBirthCounties <- getNames(subDist, "BCounty")

# Change a few manually
subBirthCounties <- subBirthCounties[-c(16, 49, 80)]

# Edit allDist to exclude the excluded isolates
newsubDist <- subDist[-c(16,49,80),-c(16,49,80)]

# Get the within and between distances for counties
subbirthCountyList <- getWithinBetween(newsubDist, subBirthCounties, FALSE)

# Plot the within between plot
plotWB(subbirthCountyList, "birth counties")

# Make VNTR colors 
subVNTRTips <- makeVNTRCols(subVNTR)

# Plot VNTR tree
plotVNTRTrees(subbedTree, subVNTRTips)

# Create an example
exampleTree <- extract.clade(subbedTree, node = 119)

# Make VNTR colours for the thing
exampleTips <- makeVNTRCols(c("2","1","2","1","1","1","1","1","1","2"))

# Get and set the margins
currentMar <- par()$mar
par(mar=c(0,0,0,0), fig=c(0,1,0,1))

# Plot tree
plot.phylo(exampleTree, edge.width = 2, font = 1, label.offset = 0.2,
           align.tip.label = FALSE, type="phylogram", cex = 0.7, show.tip.label = FALSE,
           col="grey50")

tiplabels(pch = 17, col = exampleTips,  cex = 2.5)

# Add the SNP scale
add.scale.bar(x=3,y=2,cex = 1.0, lwd = 2)
text(x=5.5,y=2.25, "SNPs")

# Reset the margins
par(mar=currentMar)



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

# Function to generate colours based on region
makeIrishRegionColours <- function(realNames){
  
  # Copy the name vector
  colourVec <- realNames
  
  # Loop thru each name
  for(index in 1:length(colourVec)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Donegal", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Louth", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Westmeath", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Cavan", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black" 
    } else if(grepl("Ireland", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if(grepl("Leitrim", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deepskyblue3"    
    } else if(grepl("Kerry", colourVec[index]) == TRUE){
        
      colourVec[index] <- "darkorange2"
    } else if(grepl("Clare", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if(grepl("Wexford", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if(grepl("Derry", colourVec[index]) == TRUE){
        
      colourVec[index] <- "black"
    } else if(grepl("Tyrone", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("Monaghan", colourVec[index]) == TRUE){
        
      colourVec[index] <- "black"
    } else if (grepl("Meath", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
      
    } else if(grepl("Sligo", colourVec[index]) == TRUE){
        
      colourVec[index] <- "deepskyblue3"
    } else if(grepl("Laois", colourVec[index]) == TRUE){
        
      colourVec[index] <- "red"
    } else if(grepl("Tipperary", colourVec[index]) == TRUE){
        
        colourVec[index] <- "darkorange2"
    } else if(grepl("Kildare", colourVec[index]) == TRUE){
          
        colourVec[index] <- "red"
    } else if(grepl("Mayo", colourVec[index]) == TRUE){
          
        colourVec[index] <- "deepskyblue3"
    } else if(grepl("Limerick", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if(grepl("Wicklow", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if(grepl("Cork", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if(grepl("Waterford", colourVec[index]) == TRUE){
      
      colourVec[index] <- "darkorange2"
    } else if (grepl("MAPK10", colourVec[index]) == TRUE){
      
      colourVec[index] <- "grey40"
    } else if (grepl("Dublin", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if (grepl("Laois", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"
    } else if (grepl("Kilkenny", colourVec[index]) == TRUE){
      
      colourVec[index] <- "red"

    } else if (grepl("Galway", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deepskyblue3"
    } else if (grepl("Roscommon", colourVec[index]) == TRUE){
      
      colourVec[index] <- "deepskyblue3"
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
    } else if(grepl("SRR", tiplabel[index]) == TRUE){
      
      dropVector <- append(dropVector, index)
    }
  }
  return(dropVector)
}

# Function to plot Irish tree
plotIrishTree <- function(tree, tipcols){
  
  # Set margins to nothing
  currentMar <- par()$mar
  par(mar=c(0,0,0,0))
  
  # Plot the national tree
  plot.phylo(tree, edge.width = 0.2, font = 1, label.offset = 0.2, 
             tip.color = tipcols,
             align.tip.label = FALSE, type="phylogram", cex = 0.5)
  
  # Add a legend
  legend("right", legend = c("Leinster", "Connaught", "Ulster", "Munster"), 
         text.col = c("red", "deepskyblue3", "black", "darkorange2"), bty = "n", cex = 1.0,
         y.intersp = 0.5)
  
  # Add the SNP scale
  add.scale.bar(x=140, y=10, cex = 2.0)
  
  # Add text to show SNPs
  text("SNPs",x=150, y=7, cex = 1.0)
  
  # Reset the margins
  par(mar=currentMar)
}

# Function to acquire polygon coordinates for map plot
getPolygonCoords <- function(counties, path) { # input is county list
  
  polygonCoords <- list() # empty list to store gps data for each county
  
  for(index in 1:length(counties)) { # for each county
    
    fileName <- paste(path, "PolygonCoords_", counties[index], ".txt", sep="" )
    # paste general path to counties, and the .txt ending to get path for each
    
    # read in the values for each county into the appropriate list segment
    polygonCoords[[counties[index]]] <- read.table(fileName, 
                                                   header = TRUE, sep = "\t")
  }
  return(polygonCoords)
}

# Function to get map limits
mapLimits <- function(polygonCoords, counties) {
  
  # Initialise vectors to store the mins and maxes of each county
  minX <- rep(NA, length(counties))
  maxX <- rep(NA, length(counties))
  minY <- rep(NA, length(counties))
  maxY <- rep(NA, length(counties))
  
  # For each county, calculate the mins and maxes, and populate the empty vectors
  for(index in 1:length(counties)) {
    
    minX[index] <- min(polygonCoords[[counties[index]]][, "X"])
    maxX[index] <- max(polygonCoords[[counties[index]]][, "X"])
    minY[index] <- min(polygonCoords[[counties[index]]][, "Y"])
    maxY[index] <- max(polygonCoords[[counties[index]]][, "Y"])
  }
  
  # Get the overall mins and maxes of X and Y to be able to set plot limits
  ranges <- c(min(minX), max(maxX), min(minY), max(maxY))
  
  return(ranges)
}

# Function to plot small map of Ireland
smallMap <- function(polygonCoords, counties, ranges) {
  
  # Create empty plot, with input of limits from above function
  plot(x=NA, y=NA,
       xlim = c(ranges[1], ranges[2]), 
       ylim = c(ranges[3], ranges[4]),
       main = "", xlab = "", ylab = "",
       bty = "n", axes = FALSE)
  
  
  for(index in 1:length(counties)) {
    
    # Check what colour needs to be assigned 
    if(counties[index] == "Donegal" ||counties[index] == "Derry" 
       ||counties[index] == "Monaghan" ||counties[index] == "Cavan"||counties[index] == "Fermanagh"
       ||counties[index] == "Tyrone"||counties[index] == "Armagh"||counties[index] == "Down"
       ||counties[index] == "Antrim"){
      
      colour = "gray20"
    } else if(counties[index] == "Mayo" || counties[index] =="Roscommon" 
              ||counties[index] == "Sligo" ||counties[index] == "Leitrim"
              ||counties[index] == "Galway"){
      
      colour = "deepskyblue3"
    } else if(counties[index] == "Clare" ||counties[index] == "Kerry" 
              ||counties[index] == "Cork" ||counties[index] == "Limerick" 
              ||counties[index] == "Tipperary" ||counties[index] == "Waterford"){
      
      colour = "darkorange2"
    } else {
      
      colour = "red"
    }
    
    # Plot the polygon of the current county
    polygon(x=polygonCoords[[counties[index]]][, "X"],
            y=polygonCoords[[counties[index]]][, "Y"],
            col = colour)
  }
}

# Function to plot Irish fan tree with corner map
plotIrishFan <- function(tree, tipcols, polygonCoords, counties, ranges){
  
  # Set margins to nothing and set figure parameters
  currentMar <- par()$mar
  par(mar=c(0,0,0,0), fig=c(0,1,0,1))
  
  # Plot national tree as a fan
  plot.phylo(tree, edge.width = 2.5, font = 1, label.offset = 0.2, 
             tip.color = tipcols, edge.color = "grey50",
             align.tip.label = FALSE, type="fan", cex = 0.5, show.tip.label = FALSE)
  
  #Add shaped tip labels
  tiplabels(pch = 18, col = tipcols,  cex = 2.3)
  
  # Add the SNP scale
  add.scale.bar(x=-100, y=-110,cex = 1.5, lcol = "grey50", lwd = 3)
  text(x=-70,y=-120, "SNPs", cex = 3)
  
  # Add a legend
  #legend(x=-110,y=100, legend = c("Leinster", "Connaught", "Ulster", "Munster"), 
  #       text.col = c("red", "deepskyblue3", "black", "darkorange2"), bty = "n", cex = 2,
  #       y.intersp = 0.6)
  
  # Set figure parameters to top right corner 
  par(fig=c(0.8,1,0.8,1), new=T)
  
  # Plot the map in top right - REQUIRES OTHER SCRIPT
  smallMap(polygonCoords, counties, ranges)
  
  # Reset the margins
  par(mar=currentMar)
}

# Function to plot unrooted tree
plotIrishUnrooted <- function(tree, tipcols, polygonCoords, counties, ranges){
  
  # Set margins to nothing and set figure parameters
  currentMar <- par()$mar
  par(mar=c(0,0,0,0), fig=c(0,1,0,1))
  
  # Plot national tree as a fan
  plot.phylo(tree, edge.width = 2.5, font = 1, label.offset = 0.2, 
             tip.color = tipcols, edge.color = "grey50",
             align.tip.label = FALSE, type="unrooted", cex = 0.5, show.tip.label = FALSE)
  
  #Add shaped tip labels
  tiplabels(pch = 18, col = tipcols,  cex = 2.3)
  
  # Add the SNP scale
  add.scale.bar(x=0, y=15, cex = 1.5, lcol = "grey50", lwd = 3)
  text(x=30,y=0, "SNPs", cex = 3)
  
  # Add a legend
  #legend(x=-110,y=100, legend = c("Leinster", "Connaught", "Ulster", "Munster"), 
  #       text.col = c("red", "deepskyblue3", "black", "darkorange2"), bty = "n", cex = 2,
  #       y.intersp = 0.6)
  
  # Set figure parameters to top right corner 
  par(fig=c(0.8,1,0.1,0.3), new=T)
  
  # Plot the map in top right - REQUIRES OTHER SCRIPT
  smallMap(polygonCoords, counties, ranges)
  
  # Reset the margins
  par(mar=currentMar)
}

# Function to plot problem herd
plotProblemHerd <- function(tree){
  
  # Get and set the margins
  currentMar <- par()$mar
  par(mar=c(0,0,0,0), fig=c(0,1,0,1))
  
  # Plot tree
  plot.phylo(tree, edge.width = 2, font = 1, label.offset = 0.2,
             align.tip.label = FALSE, type="phylogram", cex = 0.7, show.tip.label = FALSE,
             col="grey50")
  
  tiplabels(pch = 17, col = "darkorange3",  cex = 2.5)
  
  # Add the SNP scale
  add.scale.bar(x=3,y=10,cex = 1.0, lwd = 2)
  text(x=4,y=9, "SNP")
  
  # Reset the margins
  par(mar=currentMar)
}

# Function to pull out within and between SNP distances
getWithinBetween <- function(mat, name, shuffle, originalWB){
  
  if(shuffle == TRUE){
    
    shuffler <- sample(name, replace = FALSE)
    
    # Create vectors to store info in using length of original
    within <- rep(NA, length(originalWB[[1]]))
    between <- rep(NA, length(originalWB[[2]]))
    
    # Create counters for within and between
    countW <- 0
    countB <- 0
    
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
            
            countW <- countW + 1
            
            # Store distance value
            within[countW] <- mat[row,col]
            
          }else{
            
            countB <- countB + 1
            
            # Store distance value
            between[countB] <- mat[row,col]
          }
          
        }
        
      }
    }
  
  } else {
    
    shuffler <- name
    
    # Create vectors to store info in
    within <- c()
    between <- c()
    
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
            
            # Store distance value
            within <- append(within, mat[row,col])
            
          }else{
            
            # Store distance value
            between <- append(between, mat[row,col])
          }
          
        }
        
      }
    }  
  }
  return(list(w=within, b=between))        
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

# Function to plot within between plot for a list
plotWB <- function(listo, labelo){
  
  # Plot a boxplot comparing within and between
  boxplot(listo$w, listo$b, 
          main = paste("SNP distances within & between", labelo), 
          names = c("Within", "Between"),
          ylab = "SNP Difference",
          las = 1)
  stripchart(listo$w, add = TRUE, at =1, 
             method = "jitter", vertical = TRUE, col = alpha("blue",0.4),
             pch = 4)
  stripchart(listo$b, add = TRUE, at =2, 
             method = "jitter", vertical = TRUE, col = alpha("forestgreen",0.4),
             pch = 4)
  
}

# Function to run permutation and plot result
runplotPer <- function(listo, labelo, mat, name, shuffle, num, numBreaks){
  
  # Median difference of within/between groups
  diffMedWB <- median(listo$b) - median(listo$w)
  
  # Create a vector to store num values
  medVector <- rep(NA, num)
  
  # Get the within/between distances for num permuted sets
  for(run in 1:num){
    distRunnerList <- getWithinBetween(mat, name, shuffle, listo)
    
    # Get the difference of the means
    medVector[run] <- median(distRunnerList$b) - median(distRunnerList$w) 
    
  }
  
  # Plot as histogram
  xmin <- min(medVector, diffMedWB)
  xmax <- max(medVector, diffMedWB)
  
  quantiles <- quantile(medVector, c(0.025, 0.975))
  
  h <- hist(medVector, breaks=numBreaks, plot=FALSE)
  
  cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))
  
  plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
       main=paste("Clustering by", labelo), xlim=c(xmin, xmax), cex.axis=0.8, las=1)
  lines(c(diffMedWB,diffMedWB), c(0, max(h$counts)), col="blue", lwd=3)
  text(xmax-(xmax/4), num/10, cex = 0.9,
       paste("Actual Value\n= ", round(diffMedWB, digits=2)), col="blue")
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

# Function to plot VNTR tree
plotVNTRTrees <- function(tree, tipCols){
  
  # Plot the tree
  plot.phylo(tree, edge.width = 0.2, font = 1, label.offset = 0.2, 
             tip.color = tipCols,
             align.tip.label = FALSE, type="phylogram", cex = 0.2)
  
  # Add the SNP scale
  add.scale.bar(cex = 3)
  
  # Plot international tree as a fan
  plot.phylo(tree, edge.width = 2, font = 1, 
             tip.color = tipCols, edge.color = "grey50",
             align.tip.label = FALSE, type="fan", cex = 0.5, show.tip.label = FALSE)
  
  #Add shaped tip labels
  tiplabels(pch = 18, col = tipCols,  cex = 2)
  
  # Add the SNP scale
  add.scale.bar(x=30,y=-120, cex = 3,lcol = "grey50", lwd = 3)
  text(x=60,y=-130, "SNPs", cex = 3)
  
  # Add a legend
  legend(x=40, y=140, legend = c("1", "2", "3", "13", "116"), 
         text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), bty = "n", cex = 2,y.intersp = 0.6)
  
}

# Function to subsample herds
herdSubsampler <- function(herdNames, tips){
  
  # Vector to store indexes of subsampled
  indexes <- c()
  
  # Loop thru each herd
  for(herd in unique(herdNames)){
    
    # Check if the herd has 2 samples or less
    if(sum(herdNames == herd) <= 2){
      
      # Pull out the indexes and place in indexes vector
      indexes <- append(indexes, which(herdNames %in% herd))
    
    } else if(sum(herdNames == herd) >= 2){
      
      # Subsample from the big herd
      subs <- sample(which(herdNames %in% herd), 2, replace = FALSE)
      
      # Place subsample indexes in
      indexes <- append(indexes, subs)
    }
  }
  
  # Create vector of indexes from tips 
  tipper <- c(1:length(tips))
  
  # Get rid of indexes from tips
  tipsOut <- tipper[-indexes]
  
  return(tipsOut)
  
}