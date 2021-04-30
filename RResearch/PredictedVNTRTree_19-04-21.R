# Created 18-02-20
# Script to process and visualise phylogenetic trees
# Please run this FIRST, as it is a PREREQUISITE for the other scripts to work

# Load packages
library(ape)
library(phytools)
library(scales)

# Set path variables
predpathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaOct2020Format-PredictedVNTR.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/RAxML_bipartitions.RaxML-R_04-03-21"
predpathNI <- "C:/Users/UCD/Documents/Lab/CVRL MAP/NIMetaOct2020-PredictedVNTR.csv"

# Read in table of bryant isolates
isoBryantTable <- read.table(pathBryantIso,
                             header = TRUE,
                             sep = ",",
                             stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                             check.names=FALSE) # Names left as they are, no dots inserted

# Read in table of CVRL isolates
pisoCVRLTable <- read.table(predpathNewIso,
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)

# Read in table of NI isolates
pisoNITable <- read.table(predpathNI,
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

# Re-root
reRoot <- root(TheTree, node = 297)
#rotatedTree <- rotateNodes(reRoot, c(409))

# Get the Bryant names
realNames <- getBryantLabels(isoBryantTable, reRoot)

# Update the names in the tree
reRoot$tip.label <- realNames

# Get metadata for CVRL isolates
prealNames <- getCVRLLabels(pisoCVRLTable, reRoot)

# Update the names in the tree
reRoot$tip.label <- prealNames

# Get metadata for NI isolates
prealNames <- getNILabels(pisoNITable, reRoot)

# Update names
reRoot$tip.label <- prealNames

# Get rid of the non-irish isolates
dropper <- toDropInternationalTips(reRoot$tip.label)
irishOnlytree <- drop.tip(reRoot, dropper)

# Get rid of non-relevant tips
dropem <- c(17,19,20,95,96,135,136,183,184)
onlytree <- drop.tip(irishOnlytree, dropem)

# Convert branch lengths to SNP values
onlytree$edge.length <- round(onlytree$edge.length * 7843)
reRoot$edge.length <- round(reRoot$edge.length * 7843)

# Find the distances between all isolates
pallDist <- cophenetic(onlytree)
euDist <- cophenetic(reRoot)

# Round the distances
pallDist <- round(pallDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(pallDist)){
  
  pallDist[index, index] <- NA
}

# Get the herd names
herdNames <- getNames(allDist, "Herd")
euHerds <- getNames(euDist, "Herd")

# Get current county, birth county names and sameness
countyNames <- getNames(allDist, "CCounty")
birthcountyNames <- getNames(allDist, "BCounty")
sameness <- getNames(allDist, "Same")

# Get VNTR names
pvntrNames <- getNames(pallDist, "VNTR")

# Make VNTR colours
vntrTips <- makeVNTRCols(vntrNames)

# Make EU colours
euCols <- makeRegionColours(reRoot$tip.label)

# Simplify the labels
simpleLabels <- deconstructLabels(onlytree$tip.label, counties, shortCounties)

# Assign simple labels
onlytree$tip.label <- simpleLabels

#### Tree plotting (.png) ####

# Save plot as .png file (Ireland)
outputFile <- paste("VNTR_Tree_29-01-21.png", sep="")
png(outputFile, height=5000, width=4500)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 11, font = 1, label.offset = 0.2, 
           tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 2.5, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(337,308,299,284,273,264,246,212), frame = "n", cex=15, adj = c(1,1), col = "red")


# Add the SNP scale
add.scale.bar(x=10, y = 12, cex = 8, lwd = 17)
text(x=40, y =12, cex = 8, "SNPs")

# Add a legend
legend(x=9, y=210, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types")

dev.off()

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

# Save plot as .pdf file (Ireland)
outputFile <- paste("VNTR_Tree_04-03-21.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot VNTR tree
plot.phylo(onlytree, edge.width = 11, font = 1, label.offset = 0.2, tip.color = vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 2.5, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(317,285,267,259,249,222,213,199), frame = "n", cex=15, adj = c(1,1), col = "red")

# Add the SNP scale
add.scale.bar(x=10, y = 5, cex = 8, lwd = 15)
text(x=40, y =5, cex = 8, "SNPs")

# Add a legend
legend(x=9, y=160, legend = c("(42332228) - 1", "(32332228) - 2", "(32332218) - 3", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "darkorange3", "black", "darkgreen"), 
       bty = "n", cex = 10, y.intersp = 0.8, title = "INMV Types", title.col = "black")

dev.off()

# Make EU plot pdf
outputFile <- paste("EU-Tree_04-03-21.pdf", sep="")
pdf(outputFile, height=75, width=75)

# Plot EU tree
plot.phylo(reRoot, edge.width = 10, font = 1, label.offset = 0.2, 
           show.tip.label = FALSE,
           align.tip.label = FALSE, type="phylogram", cex = 30, no.margin = TRUE)
nodelabels(text= c("A","B","C","D","E","F","G","H"), node = c(472,425,403,391,379,334,303,296), frame = "n", cex=15, adj = c(1,0), col = "red")


#Add shaped tip labels
tiplabels(pch = 18, col = euCols,  cex = 9)

# Add the SNP scale
add.scale.bar(x=100, y=0, cex = 10, lwd = 15)
text(x=167, y=0, cex = 10, "SNPs")

# Add a legend
legend(x=150, y=280, legend = c("Ireland", "UK", "England", "Scotland", "Wales",
                                "Italy", "Spain", "France", "Germany", "Netherlands",
                                "Czech Rep.", "Greece", "Norway"), 
       text.col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
                    "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
                    "mediumblue", "slateblue", "purple"), pch = 18,
       col = c("darkgreen", "firebrick4", "lightpink2", "steelblue3", "deeppink",
               "aquamarine2", "goldenrod3", "royalblue4", "black", "orangered",
               "mediumblue", "slateblue", "purple"),
       bty = "n", cex = 10, y.intersp = 0.8, title = "Country", title.col = "black")

dev.off()

#### Code for clustering analyses of herds/vntrs/counties ####

pvdistList <- getWithinBetween(pallDist, pvntrNames, FALSE)

# Get the difference of the means
pdiffWB <- median(pvdistList$Between) - median(pvdistList$Within) 

# Plot a boxplot comparing within and between
boxplot(pvdistList$Within, pvdistList$Between, 
        main = "SNP distances within & between predicted VNTR types", 
        names = c("Within", "Between"),
        ylab = "SNP Difference")
stripchart(pvdistList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("lightblue",0.4),
           pch = 4)
stripchart(pvdistList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("lightgreen",0.4),
           pch = 4)



# Create an empty vector
medianVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  distRunnerList <- getWithinBetween(allDist, vntrNames, TRUE)
  
  # Get the difference of the means
  medianVector[run] <- median(distRunnerList$Between) - median(distRunnerList$Within) 
  
}

# Plot as histogram
xmin <- min(medianVector, diffWB)
xmax <- max(medianVector, diffWB)

quantiles <- quantile(medianVector, c(0.025, 0.975))

h <- hist(medianVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Distribution of differences between medians", xlim=c(xmin, 25), cex.axis=0.8, las=1)
lines(c(diffWB,diffWB), c(0, max(h$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffWB, digits=2)), col="blue")

cdistList <- getWithinBetween(allDist, countyNames, FALSE)

# Get the difference of the means
cdiffWB <- median(cdistList$Between) - median(cdistList$Within) 

# Plot a boxplot comparing within and between
boxplot(cdistList$Within, cdistList$Between, 
        main = "SNP distances within & between counties", 
        names = c("Within", "Between"),
        ylab = "SNP Difference")
stripchart(cdistList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("lightblue",0.4),
           pch = 4)
stripchart(cdistList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("lightgreen",0.4),
           pch = 4)



# Create an empty vector
cmedianVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  cdistRunnerList <- getWithinBetween(allDist, countyNames, TRUE)
  
  # Get the difference of the means
  cmedianVector[run] <- median(cdistRunnerList$Between) - median(cdistRunnerList$Within) 
  
}

# Plot as histogram
cxmin <- min(cmedianVector, cdiffWB)
cxmax <- max(cmedianVector, cdiffWB)

cquantiles <- quantile(cmedianVector, c(0.025, 0.975))

ch <- hist(cmedianVector, breaks=30, plot=FALSE)

ccuts <- cut(ch$breaks, c(-Inf, cquantiles[1], cquantiles[2], Inf))

plot(h, col=c("red", "white", "red")[ccuts], xlab="Difference",
     main="Distribution of differences between medians", xlim=c(cxmin, 25), cex.axis=0.8, las=1)
lines(c(cdiffWB,cdiffWB), c(0, max(ch$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(cdiffWB, digits=2)), col="blue")

hdistList <- getWithinBetween(allDist, herdNames, FALSE)

# Get the difference of the means
hdiffWB <- median(hdistList$Between) - median(hdistList$Within) 

# Plot a boxplot comparing within and between
boxplot(hdistList$Within, hdistList$Between, 
        main = "SNP distances within & between herds", 
        names = c("Within", "Between"),
        ylab = "SNP Difference")
stripchart(hdistList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("lightblue",0.4),
           pch = 4)
stripchart(hdistList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("lightgreen",0.4),
           pch = 4)



# Create an empty vector
hmedianVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  hdistRunnerList <- getWithinBetween(allDist, herdNames, TRUE)
  
  # Get the difference of the means
  hmedianVector[run] <- median(hdistRunnerList$Between) - median(hdistRunnerList$Within) 
  
}

# Plot as histogram
hxmin <- min(hmedianVector, hdiffWB)
hxmax <- max(hmedianVector, hdiffWB)

hquantiles <- quantile(hmedianVector, c(0.025, 0.975))

hh <- hist(hmedianVector, breaks=30, plot=FALSE)

hcuts <- cut(hh$breaks, c(-Inf, hquantiles[1], hquantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Distribution of differences between medians", xlim=c(hxmin, 25), cex.axis=0.8, las=1)
lines(c(hdiffWB,hdiffWB), c(0, max(hh$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(hdiffWB, digits=2)), col="blue")



# Get the difference of the means
meandiffWB <- mean(distList$Between) - mean(distList$Within) 

# Create an empty vector
meanVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  meandistRunnerList <- getWithinBetween(allDist, vntrNames, TRUE)
  
  # Get the difference of the means
  meanVector[run] <- mean(meandistRunnerList$Between) - mean(distRunnerList$Within) 
  
}

# Plot as histogram
xmin <- min(meanVector, meandiffWB)
xmax <- max(meanVector, meandiffWB)

mquantiles <- quantile(meanVector, c(0.025, 0.975))

mh <- hist(meanVector, breaks=30, plot=FALSE)

mcuts <- cut(mh$breaks, c(-Inf, mquantiles[1], mquantiles[2], Inf))

plot(mh, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Distribution of differences between means", xlim=c(xmin, 25), cex.axis=0.8, las=1)
lines(c(meandiffWB,meandiffWB), c(0, max(h$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(meandiffWB, digits=2)), col="blue")

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

