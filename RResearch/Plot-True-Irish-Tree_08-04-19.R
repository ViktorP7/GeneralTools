# 11/03/19 modified for a RAxML tree input file
# 25-03-19 modified for extra isolates
# 08-04-19 modified for better root and to find VNTR group distances
# 07-05-19 modified for new metadata
# 12-06-19 modified for new data

# Load packages
library(ape)
library(phytools)
library(scales)

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

pdf("June2019Tree.pdf", width=20, height=20)

# Root the tree at 466 - ancestor rooted in the Bryant paper
rootree <- root(TheTree, node = 466)

# Drop tips for far away ancestors (silvaticum and hominissius)
dropNumbers <- c(186,187)
droppedTree <- drop.tip(rootree, dropNumbers)

# Extract the clade that doesn't have all the distant sheep and cows
extractedTree <- extract.clade(droppedTree, node = 452)

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(extractedTree$tip.label)
irishOnlytree <- drop.tip(extractedTree, dropper)

# Convert branch lengths to SNP values
irishOnlytree$edge.length <- irishOnlytree$edge.length * 48283

# Get the rounded values o the lengths
roundedSNPs <- round(irishOnlytree$edge.length)

# Assign floored SNPs
irishOnlytree$edge.length <- roundedSNPs

# Get the colours
tipColours <- makeIrishRegionColours(irishOnlytree$tip.label)

# Plot Irish tree
plotIrishTree(irishOnlytree, tipColours)

# Plot Irish fan tree
plotIrishFan(irishOnlytree, tipColours, polygonCoords, counties)

# Extract Cork 10 herd only
cork10 <- extract.clade(irishOnlytree, node = 165)

# Get the colours
corktipColours <- makeIrishRegionColours(cork10$tip.label)

# Refresh par
par(mar=c(0,0,0,0), fig=c(0,1,0,1))

# Plot Cork 10
plot.phylo(cork10, edge.width = 0.2, font = 1, label.offset = 0.2, 
           tip.color = corktipColours,
           align.tip.label = FALSE, type="phylogram", cex = 0.7)


# Add the SNP scale
add.scale.bar(x=6,y=10,cex = 1.0)
text(x=6.5,y=9, "SNP")

dev.off()

#### Post Tree Analysis ####

# Find INMV type 1 strains
inmv1 <- findVNTRs(irishOnlytree$tip.label, "1")

# Find INMV type 2 strains
inmv2 <- findVNTRs(irishOnlytree$tip.label, "2")

# Find the distances between all isolates
allDist <- cophenetic(irishOnlytree)

# Round the distances
allDist <- round(allDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(allDist)){
  
  allDist[index, index] <- NA
}

# other INMV locations
others <- c(3,69,121)

# Combine INMV vectors and sort
allINMV <- c(inmv1, inmv2, others)
sortedAll <- sort(allINMV)

# Filter out all INMVs
inmvDist <- allDist[allINMV, allINMV]

#### No need for these below ####

# Filter out type 1 INMVs
inmv1Dist <- allDist[inmv1, inmv1]

# Filter out type 2 INMVs
inmv2Dist <- allDist[inmv2, inmv2]

# Find the max distances
maxINMV1 <- max(inmv1Dist, na.rm = TRUE)
maxINMV2 <- max(inmv2Dist, na.rm = TRUE)
maxAll <- max(inmvDist, na.rm = TRUE)

# Find the min distances
minINMV1 <- min(inmv1Dist, na.rm = TRUE)
minINMV2 <- min(inmv2Dist, na.rm = TRUE)
minAll <- min(inmvDist, na.rm = TRUE)

# Find the median value
medINMV1 <- median(inmv1Dist, na.rm = TRUE)
medINMV2 <- median(inmv2Dist, na.rm = TRUE)
medAll <- median(inmvDist, na.rm = TRUE)

# Find the mean value
meanINMV1 <- floor(mean(inmv1Dist, na.rm = TRUE))
meanINMV2 <- floor(mean(inmv2Dist, na.rm = TRUE))
meanAll <- floor(mean(inmvDist, na.rm = TRUE))

# Find the frequency of SNP distances
allFreqs <- findSNPFreqs(inmvDist)

# Add up all the freqs from previous rows
addedAllFreqs <- allFreqs

for(row in 2:nrow(allFreqs)){
  
  addedAllFreqs[row,1] <- allFreqs[row,1] + allFreqs[row-1,1]
  addedAllFreqs[row,2] <- allFreqs[row,2] + allFreqs[row-1,2]
  addedAllFreqs[row,3] <- allFreqs[row,3] + allFreqs[row-1,3]
}

# Get proportions
propAll <- processFreqs(allFreqs)
propAdded <- processFreqs(addedAllFreqs)

# Plot the proportions - don't actually use this
plot(propAll[,1], col = "blue", ylim = c(0,100),
     ylab = "Frequency", xlab = "SNPs Difference", xaxt = "n")
axis(1, at=1:10, labels = rownames(propAll), las = 2)
points(propAll[,2], col = "red")
lines(propAll[,1], col = "red")
lines(propAll[,2], col = "blue", lty = 2)

#### Continue using these below for VNTR ####

# Process distance matrices for each VNTR type
inmvDist[upper.tri(inmvDist)] <- NA

# Get names of matrix
nameVec <- getNames(inmvDist)

# Get the within and between distances
distList <- getWithinBetween(inmvDist, nameVec, FALSE)

# Plot a boxplot comparing within and between
boxplot(distList$Within, distList$Between, 
        main = "SNP distances within & between VNTR types", 
        names = c("Within", "Between"),
        ylab = "SNP Difference",
        las = 1)
stripchart(distList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",0.4),
           pch = 4)
stripchart(distList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("forestgreen",0.4),
           pch = 4)

# Do with median differences
diffMedWB <- median(distList$Between) - median(distList$Within)

# Create a vector to store 10k values
medVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
system.time(for(run in 1:10000){
  distRunnerList <- getWithinBetween(inmvDist, nameVec, TRUE)
  
  # Get the difference of the means
  medVector[run] <- median(distRunnerList$Between) - median(distRunnerList$Within) 
  
})

# Plot as histogram
xmin <- min(medVector, diffMedWB)
xmax <- max(medVector, diffMedWB)

quantiles <- quantile(medVector, c(0.025, 0.975))

h <- hist(medVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Isolate VNTR Type", xlim=c(xmin, 10), cex.axis=0.8, las=1)
lines(c(diffMedWB,diffMedWB), c(0, max(h$counts)), col="blue", lwd=3)
text(8, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffMedWB, digits=2)), col="blue")

#### Do with herd names now instead of VNTR, get herds first ####
allHerds <- isolateHerder(irishOnlytree$tip.label)

# Get herd distances
herdDist <- allDist[allHerds, allHerds]

# Process distance matrices for each herd
herdDist[upper.tri(herdDist)] <- NA

# Get the herd names
herdNames <- getHerdNames(herdDist)

# Get the within and between distances for herds
distHerdList <- getWithinBetween(herdDist, herdNames, FALSE)

# Get the difference of the medians
diffMedHerdWB <- median(distHerdList$Between) - median(distHerdList$Within)

# Plot a boxplot comparing within and between for herds
boxplot(distHerdList$Within, distHerdList$Between, 
        main = "SNP distances within & between herds", 
        names = c("Within", "Between"),
        ylab = "SNP Difference",
        las = 1)
stripchart(distHerdList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",1),
           pch = 4)
stripchart(distHerdList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("forestgreen",0.4),
           pch = 4)

# Create a vector to store 10k values
medHerdVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
system.time(for(run in 1:10000){
  distHerdRunnerList <- getWithinBetween(herdDist, herdNames, TRUE)
  
  # Get the difference of the medians
  medHerdVector[run] <- median(distHerdRunnerList$Between) - median(distHerdRunnerList$Within) 
  
})

# Plot as histogram
xmin <- min(medHerdVector, diffMedHerdWB)
xmax <- max(medHerdVector, diffMedHerdWB)

quantiles <- quantile(medHerdVector, c(0.025, 0.975))

h <- hist(medHerdVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Isolate Herd", xlim=c(xmin, 190), cex.axis=0.8, las=1)
lines(c(diffMedHerdWB,diffMedHerdWB), c(0, max(h$counts)), col="blue", lwd=3)
text(165, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffMedHerdWB, digits=2)), col="blue")

# Now do the same but look at the county level instead of herd
countyNames <- processHerds(herdNames)
# Change a few manually
countyNames[c(84,85,118, 119, 120, 74)] <- "Cork"
countyNames[c(8,15)] <- "North1"
countyNames[c(142,143)] <- "North2"

# Get the within and between distances for counties
distCountyList <- getWithinBetween(herdDist, countyNames, FALSE)

# Get the difference of the medians
diffMedCountyWB <- median(distCountyList$Between) - median(distCountyList$Within)

# Plot a boxplot comparing within and between for counties
boxplot(distCountyList$Within, distCountyList$Between, 
        main = "SNP distances within & between counties", 
        names = c("Within", "Between"),
        ylab = "SNP Difference",
        las = 1)
stripchart(distCountyList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",0.6),
           pch = 4)
stripchart(distCountyList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("forestgreen",0.4),
           pch = 4)

# Create a vector to store 10k values
medCountyVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
system.time(for(run in 1:10000){
  distCountyRunnerList <- getWithinBetween(herdDist, countyNames, TRUE)
  
  # Get the difference of the medians
  medCountyVector[run] <- median(distCountyRunnerList$Between) - median(distCountyRunnerList$Within) 
  
})

# Plot as histogram
xmin <- min(medCountyVector, diffMedCountyWB)
xmax <- max(medCountyVector, diffMedCountyWB)

quantiles <- quantile(medCountyVector, c(0.025, 0.975))

h <- hist(medCountyVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Isolate County", xlim=c(xmin, 80), cex.axis=0.8, las=1)
lines(c(diffMedCountyWB,diffMedCountyWB), c(0, max(h$counts)), col="blue", lwd=3)
text(65, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffMedCountyWB, digits=2)), col="blue")

# Look at the birth county now
birthcountyNames <- getBirthCounties(herdDist)
# Change a few manually
birthcountyNames <- birthcountyNames[-c(8,9,15,52,74,84,85,90,117,118,119,120,142,143)]

newHerdDist <- herdDist[-c(8,9,15,52,74,84,85,90,117,118,119,120,142,143),-c(8,9,15,52,74,84,85,90,117,118,119,120,142,143)]

# Get the within and between distances for counties
distbirthCountyList <- getWithinBetween(newHerdDist, birthcountyNames, FALSE)

# Get the difference of the medians
diffMedbirthCountyWB <- median(distbirthCountyList$Between) - median(distbirthCountyList$Within)

# Plot a boxplot comparing within and between for counties
boxplot(distbirthCountyList$Within, distbirthCountyList$Between, 
        main = "SNP distances within & between birth counties", 
        names = c("Within", "Between"),
        ylab = "SNP Difference",
        las = 1)
stripchart(distbirthCountyList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",0.6),
           pch = 4)
stripchart(distbirthCountyList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("forestgreen",0.4),
           pch = 4)

# Create a vector to store 10k values
medbirthCountyVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
system.time(for(run in 1:10000){
  distbirthCountyRunnerList <- getWithinBetween(newHerdDist, birthcountyNames, TRUE)
  
  # Get the difference of the medians
  medbirthCountyVector[run] <- median(distbirthCountyRunnerList$Between) - median(distbirthCountyRunnerList$Within) 
  
})

# Plot as histogram
xmin <- min(medbirthCountyVector, diffMedbirthCountyWB)
xmax <- max(medCountyVector, diffMedCountyWB)

quantiles <- quantile(medbirthCountyVector, c(0.025, 0.975))

h <- hist(medbirthCountyVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Isolate Birth County", xlim=c(xmin, 50), cex.axis=0.8, las=1)
lines(c(diffMedbirthCountyWB,diffMedbirthCountyWB), c(0, max(h$counts)), col="blue", lwd=3)
text(35, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffMedbirthCountyWB, digits=2)), col="blue")
###########################################################################
###FUNCTIONS###FUNCTIONS###FUNCTIONS###FUNCTIONS###FUNCTIONS###FUNCTIONS###
###########################################################################

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
    if(grepl("NIreland", colourVec[index]) == TRUE){
      
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
    } else if(grepl("NA", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"

    } else if(grepl("Donegal", colourVec[index]) == TRUE){
      
      colourVec[index] <- "black"
    } else if(grepl("CIT", colourVec[index]) == TRUE){
      
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
    } 
  }
  return(dropVector)
}

# Function to plot Irish tree
plotIrishTree <- function(tree, tipcols){
  
  # Set margins to nothing
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
  add.scale.bar(cex = 2.0)
}

# Function to plot Irish fan tree with corner map
plotIrishFan <- function(tree, tipcols, polygonCoords, counties){
  
  # Set margins to nothing and set figure parameters
  par(mar=c(0,0,0,0), fig=c(0,1,0,1))
  
  # Plot national tree as a fan
  plot.phylo(tree, edge.width = 2.5, font = 1, label.offset = 0.2, 
             tip.color = tipcols, edge.color = "grey50",
             align.tip.label = FALSE, type="fan", cex = 0.5, show.tip.label = FALSE)
  
  #Add shaped tip labels
  tiplabels(pch = 18, col = tipcols,  cex = 2.5)
  
  # Add the SNP scale
  add.scale.bar(cex = 1.5, lcol = "grey50", lwd = 3)
  text(x=-70,y=-120, "SNPs", cex = 3)
  
  # Add a legend
  legend(x=-110,y=100, legend = c("Leinster", "Connaught", "Ulster", "Munster"), 
         text.col = c("red", "deepskyblue3", "black", "darkorange2"), bty = "n", cex = 2,
         y.intersp = 0.6)
  
  # Set figure parameters to top right corner 
  par(fig=c(0.8,1,0.8,1), new=T)
  
  # Plot the map in top right - REQUIRES OTHER SCRIPT
  smallMap(polygonCoords, counties)
}

# Function to pull out tiplabels corresponding to a VNTR type
findVNTRs <- function(tiplabel, VNTR){
  
  # Create vector of tips to find
  finderVec <- c()
  
  # Loop thru the tips
  for(index in 1:length(tiplabel)){
    
    # Split the string and take the 3rd value
    vntrInfo <- strsplit(tiplabel[index], split = "_")[[1]][4]
    
    # Skip if it's an NA
    if(is.na(vntrInfo) == TRUE){
      
      next
    
    }else if(vntrInfo == VNTR){
      
      finderVec <- append(finderVec, index)
      
    }
  }
  
  return(finderVec)
}

# Function to count frequency of SNP difference 
findSNPFreqs <- function(mat){
  
  # Initialise a matrix to store freq counts
  counter <- matrix(data = 0, nrow = 10, ncol = 3, 
                    dimnames = list(c("<5", "<10", "<20", 
                                      "<30", "<40", "<50",
                                      "<75", "<100", "<200", "<250"),
                                    c("Same", "NotSame", "Total")))
  
  # Loop thru each row of the input and count how many times each appears
  for(row in 1:nrow(mat)){
    
    vRow <- strsplit(rownames(mat)[row], split = "_")[[1]][3]
    
    for(col in 1:ncol(mat)){
      vCol <- strsplit(colnames(mat)[col], split = "_")[[1]][3]
      
      if(vRow == vCol){
        
        if(is.na(mat[row,col]) == TRUE){
          
          next
          
        }else if(mat[row,col] < 5){
          counter[1,1] <- counter[1,1] + 1
        }else if(mat[row,col] < 10){
          counter[2,1] <- counter[2,1] + 1
        }else if(mat[row,col] < 20){  
          counter[3,1] <- counter[3,1] + 1
        }else if(mat[row,col] < 30){
          counter[4,1] <- counter[4,1] + 1
        }else if(mat[row,col] < 40){
          counter[5,1] <- counter[5,1] + 1
        }else if(mat[row,col] < 50){
          counter[6,1] <- counter[6,1] + 1
        }else if(mat[row,col] < 75){
          counter[7,1] <- counter[7,1] + 1
        }else if(mat[row,col] < 100){
          counter[8,1] <- counter[8,1] + 1
        }else if(mat[row,col] < 200){
          counter[9,1] <- counter[9,1] + 1
        }else if(mat[row,col] < 250){
          counter[10,1] <- counter[10,1] + 1
        }  
      } else{
        
        if(mat[row,col] < 5){
          counter[1,2] <- counter[1,2] + 1
        }else if(mat[row,col] < 10){
          counter[2,2] <- counter[2,2] + 1
        }else if(mat[row,col] < 20){  
          counter[3,2] <- counter[3,2] + 1
        }else if(mat[row,col] < 30){
          counter[4,2] <- counter[4,2] + 1
        }else if(mat[row,col] < 40){
          counter[5,2] <- counter[5,2] + 1
        }else if(mat[row,col] < 50){
          counter[6,2] <- counter[6,2] + 1
        }else if(mat[row,col] < 75){
          counter[7,2] <- counter[7,2] + 1
        }else if(mat[row,col] < 100){
          counter[8,2] <- counter[8,2] + 1
        }else if(mat[row,col] < 200){
          counter[9,2] <- counter[9,2] + 1
        }else if(mat[row,col] < 250){
          counter[10,2] <- counter[10,2] + 1
        }  
        
      }
    }
  }
  
  # Get the totals
  for(row in 1:nrow(counter)){
    
    counter[row,3] <- counter[row,1] + counter[row,2]
  }
  
  return(counter)
}

# Function to process SNP distance matrix and give percentage proportions
processFreqs <- function(freqs){
  
  # Duplicate the input
  propFreqs <- freqs
  
  # Loop thru and get proportions
  for(row in 1:nrow(freqs)){
    
    propFreqs[row,1] <- freqs[row,1]/freqs[row,3]*100
    propFreqs[row,2] <- freqs[row,2]/freqs[row,3]*100
    propFreqs[row,3] <- 100
  }
  
  return(propFreqs)
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

# Function to pull out matrix names and simplify to VNTR types
getNames <- function(mat){
  
  # Create a pair of vectors for row and col names
  rowcolNames <- rep(NA, length(colnames(mat)))
  
  # Loop thru and split off what's needed and fill into vectors
  for(index in 1:length(rowcolNames)){
    
    # Store row name
    vRow <- strsplit(rownames(mat)[index], split = "_")[[1]][4]
    
    rowcolNames[index] <- vRow
  }
  return(rowcolNames)
}

# Function to pull out matrix names and simplify to herd names
getHerdNames <- function(mat){
  
  # Create a pair of vectors for row and col names
  rowcolNames <- rep(NA, length(colnames(mat)))
  
  # Loop thru and split off what's needed and fill into vectors
  for(index in 1:length(rowcolNames)){
    
    # Store row name
    one <- strsplit(colnames(mat)[index], split = "_")[[1]][2]
    two <- strsplit(colnames(mat)[index], split = "_")[[1]][3]
    
    vRow <- paste(one,"_",two)
    
    rowcolNames[index] <- vRow
  }
  return(rowcolNames)
}

# Function to process herd names into just county names
processHerds <- function(herds){
  
  # Create vector to store counties
  counties <- rep(NA, length(herds))
  
  # Loop thru the herdnames
  for(index in 1:length(herds)){
    
    counties[index] <- strsplit(herds[index], split = "_")[[1]][1]
  }
  return(counties)
}

# Function to pull out tiplabels corresponding to herds
isolateHerder <- function(tiplabel){
  
  # Create vector of tips to find
  finderVec <- c()
  
  # Loop thru the tips
  for(index in 1:length(tiplabel)){
    
    # Split the string and take the 3rd value
    one <- strsplit(tiplabel[index], split = "_")[[1]][2]
    two <- strsplit(tiplabel[index], split = "_")[[1]][3]
    
    herdInfo <- paste(one, "_", two)
    
    
    # Skip if it's an NA
    if(is.na(herdInfo) == TRUE){
      
      next
      
    }else{
      
      finderVec <- append(finderVec, index)
      
    }
  }
  
  return(finderVec)
}

# Function to pull out birth counties from matrix names
getBirthCounties <- function(mat){
  
  # Create a pair of vectors for row and col names
  rowcolNames <- rep(NA, length(colnames(mat)))
  
  # Loop thru and split off what's needed and fill into vectors
  for(index in 1:length(rowcolNames)){
    
    # Store row name
    vRow <- strsplit(colnames(mat)[index], split = "_")[[1]][5]
    
    rowcolNames[index] <- vRow
  }
  return(rowcolNames)
}