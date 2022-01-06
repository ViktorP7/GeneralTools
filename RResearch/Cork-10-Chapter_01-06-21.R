# Created 01-06-21

# Load packages
library(ape)
library(phytools)
library(scales)
library(adegenet)

# Set path variables
pathNewIso <- "C:/Users/UCD/Documents/Lab/CVRL MAP/MetaMay2021Format.csv"
pathBryantIso <- "C:/Users/UCD/Documents/Papers/Bryant 2016 Table S1.csv"
pathC10Tree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Cork10Fastqs/C10genie/RAxML_bipartitions.RaxML-R_01-06-21"
pathTyrTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SnipGenieTyr/RAxML_bipartitions.Tyr"
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
c10Tree <- read.tree(pathC10Tree)
tyrTree <- read.tree(pathTyrTree)

# Re-root
c10reRoot <- root(c10Tree, node = 105)
tyrReroot <- root(tyrTree, node = 48)

# Get the Bryant names
c10realNames <- altinputgetBryantLabels(isoBryantTable, c10reRoot)

# Update the names in the tree
c10reRoot$tip.label <- c10realNames

# Get metadata for CVRL isolates
c10realNames <- altgetCVRLLabels(isoCVRLTable, c10reRoot)
tyrReals <- altgetCVRLLabels(isoCVRLTable, tyrReroot)

# Update the names in the tree
c10reRoot$tip.label <- c10realNames
tyrReroot$tip.label <- tyrReals

# Get NI meta
tyrReals <- getNILabels(isoNITable, tyrReroot)
tyrReroot$tip.label <- tyrReals

#Get rid of non-relevant tips
c10reRoot <- drop.tip(c10reRoot, tip = 5)
tyrReroot <- drop.tip(tyrReroot, tip = 2)

# Convert branch lengths to SNP values
c10reRoot$edge.length <- round(c10reRoot$edge.length * 135)
tyrReroot$edge.length <- round(tyrReroot$edge.length * 103)

# Find the distances between all isolates
#allDist <- cophenetic(onlytree)
c10Dist <- cophenetic(c10reRoot)
tyrDist <- cophenetic(tyrReroot)

# Round the distances
c10Dist <- round(c10Dist)
tyrDist <- round(tyrDist)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(c10Dist)){
  
  c10Dist[index, index] <- NA
}

for(index in 1:nrow(tyrDist)){
  
  tyrDist[index, index] <- NA
}

# Get the herd names
c10herdNames <- getNames(c10Dist, "Herd")
tyrherdNames <- getNames(tyrDist, "Herd")

# Get current county, birth county names and sameness
c10countyNames <- getNames(c10Dist, "CCounty")
c10birthcountyNames <- getNames(c10Dist, "BCounty")
c10sameness <- getNames(c10Dist, "Same")

tyrcountyNames <- getNames(tyrDist, "CCounty")
tyrbirthcountyNames <- getNames(tyrDist, "BCounty")
tyrsameness <- getNames(tyrDist, "Same")


# Get VNTR names
c10vntrNames <- getNames(c10Dist, "VNTR")
tyrvntrNames <- getNames(tyrDist, "VNTR")

# Make VNTR colours
c10vntrTips <- makeVNTRCols(c10vntrNames)
tyrvntrTips <- makeVNTRCols(tyrvntrNames)


# Simplify the labels
c10simpleMat <- altdeconstructLabels(c10reRoot$tip.label, counties, shortCounties)
tyrsimpleMat <- altdeconstructLabels(tyrReroot$tip.label, counties, shortCounties)


# Assign simple labels
c10reRoot$tip.label <- c10simpleMat[1,]
tyrReroot$tip.label <- tyrsimpleMat[1,]

#### Tree plotting (.pdf) ####

# Save plot as .pdf file (Ireland)
outputFile <- paste("C10_Tree_07-10-21.pdf", sep="")
pdf(outputFile, height=50, width=50)

# Plot VNTR tree
plot.phylo(c10reRoot, edge.width = 17, font = 2, label.offset = 0.5, tip.color = c10vntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 3.5, no.margin = TRUE)

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

outputFile <- paste("Tyr_Tree_07-10-21.pdf", sep="")
pdf(outputFile, height=50, width=50)

# Plot VNTR tree
plot.phylo(tyrReroot, edge.width = 17, font = 2, label.offset = 0.5, tip.color = tyrvntrTips,
           align.tip.label = FALSE, type="phylogram", cex = 3.5, no.margin = TRUE)

tiplabels(pch = 18, frame = "n", col = tyrsimpleMat[2,], cex=10)
# Add the SNP scale
add.scale.bar(x=11, y = 1, cex = 7, lwd = 14)
text(x=17.5,y=1, cex = 7, "SNPs")



# Add a legend
legend(x=14, y=35, legend = c("(42332228) - 1", "(32332228) - 2", "(22332228) - 13", "(41332228) - 116"), 
       text.col = c("red", "deepskyblue3", "black", "darkgreen"), 
       bty = "n", cex = 9, y.intersp = 0.8, title = "INMV Types", title.col = "black")
legend(x=15, y=25, legend = c("2013","2014", "2016", "2017", "2019"), 
       text.col = c("red", "goldenrod3", "steelblue3", "palegreen3", "orchid"),
       pch = c(23,23,23,23,23), pt.bg = c("red","goldenrod3",  "steelblue3", "palegreen3", "orchid"), pt.cex = 13,
       bty = "n", cex = 9, y.intersp = 0.8, title = "Isolation Year", title.col = "black")

dev.off()

#### Check VNTR clustering ####

# Assign scottish VNTRs (both are 1)
c10vntrNames[1] <- 1
c10vntrNames[99] <- 1

distList <- getWithinBetween(c10Dist, c10vntrNames, FALSE)
tyrdistList <- getWithinBetween(tyrDist,tyrvntrNames, FALSE)

# Get the difference of the means
diffWB <- median(distList$Between) - median(distList$Within)
tyrdiffWB <- median(tyrdistList$Between) - median(tyrdistList$Within)

# Plot a boxplot comparing within and between
boxplot(distList$Within, distList$Between, tyrdistList$Within, tyrdistList$Between, 
        main = "SNP distances within & between VNTR types", 
        names = c("Within C10", "Between C10", "Within Tyr", "Between Tyr"),
        ylab = "SNP Difference")
stripchart(distList$Within, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("lightblue",0.4),
           pch = 4)
stripchart(distList$Between, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("lightgreen",0.4),
           pch = 4)
stripchart(tyrdistList$Within, add = TRUE, at =3, 
           method = "jitter", vertical = TRUE, col = alpha("lightblue",0.4),
           pch = 4)
stripchart(tyrdistList$Between, add = TRUE, at =4, 
           method = "jitter", vertical = TRUE, col = alpha("lightgreen",0.4),
           pch = 4)


range(c10Dist[which(c10vntrNames == 1),which(c10vntrNames == 1)], na.rm=T)
range(c10Dist[which(c10vntrNames == 2),which(c10vntrNames == 2)], na.rm=T)
median(c10Dist[which(c10vntrNames == 1),which(c10vntrNames == 1)], na.rm=T)
median(c10Dist[which(c10vntrNames == 2),which(c10vntrNames == 2)], na.rm=T)
range(tyrDist[which(tyrvntrNames == 1),which(tyrvntrNames == 1)], na.rm=T)
range(tyrDist[which(tyrvntrNames == 2),which(tyrvntrNames == 2)], na.rm=T)
median(tyrDist[which(tyrvntrNames == 1),which(tyrvntrNames == 1)], na.rm=T)
median(tyrDist[which(tyrvntrNames == 2),which(tyrvntrNames == 2)], na.rm=T)

# Create an empty vector
medianVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  distRunnerList <- getWithinBetween(c10Dist, c10vntrNames, TRUE)
  
  # Get the difference of the means
  medianVector[run] <- mean(distRunnerList$Between) - mean(distRunnerList$Within) 
  
}

# Plot as histogram
xmin <- min(medianVector, diffWB)
xmax <- max(medianVector, diffWB)

quantiles <- quantile(medianVector, c(0.025, 0.975))

h <- hist(medianVector, breaks=30, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf, quantiles[1], quantiles[2], Inf))

plot(h, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Distribution of differences between means C10", xlim=c(xmin, 25), cex.axis=0.8, las=1)
lines(c(diffWB,diffWB), c(0, max(h$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(diffWB, digits=2)), col="blue")

tmedianVector <- rep(NA, 10000)

# Do the same except this time shuffle 10k times
for(run in 1:10000){
  tdistRunnerList <- getWithinBetween(tyrDist, tyrvntrNames, TRUE)
  
  # Get the difference of the means
  tmedianVector[run] <- mean(tdistRunnerList$Between) - mean(tdistRunnerList$Within) 
  
}

# Plot as histogram
xmin <- min(tmedianVector, tyrdiffWB)
xmax <- max(tmedianVector, tyrdiffWB)

tquantiles <- quantile(tmedianVector, c(0.025, 0.975))

th <- hist(tmedianVector, breaks=30, plot=FALSE)

tcuts <- cut(th$breaks, c(-Inf, tquantiles[1], tquantiles[2], Inf))

plot(th, col=c("red", "white", "red")[cuts], xlab="Difference",
     main="Distribution of differences between means Tyr CaA", xlim=c(xmin, 25), cex.axis=0.8, las=1)
lines(c(tyrdiffWB,tyrdiffWB), c(0, max(th$counts)), col="blue", lwd=3)
text(20, 750, cex = 0.9,
     paste("Actual Value\n= ", round(tyrdiffWB, digits=2)), col="blue")

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

#Function to get labels for NI isolates
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