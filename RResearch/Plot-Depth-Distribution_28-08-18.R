# Load libraries/packages
library(xlsx)

# Set path variable
pathVcf <- "path/to/vcf"
pathAnno <- "path/to/annotation"

# Read in vcf
vcfPTable <- read.table(pathVcf,
                       header = FALSE,
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted
# Read in gff
gffTable <- read.table(pathAnno,
                       header = FALSE,
                       quote = "", # Avoids In scan...EOF within quoted string error
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted

# Remove unecessary row and columns from gff
gffTable <- gffTable[, -c(1,2,6,7,8)]
gffTable <- gffTable[-which(gffTable$V3 == "gene"),]
gffTable <- gffTable[-1,]
gffTable <- gffTable[, -1]

# Create vectors to store positions and quality info
positionVectorP <- vcfPTable[,"V2"]
qualityVectorP <- vcfPTable[,"V8"]

# Vector to store DP only
dpVectorP <- rep(NA, length(qualityVectorP))

# Loop thru quality vector and get rid of irrelevant information, store only DP value
for(index in 1:length(qualityVectorP)){
  
  # Create variable to store split position
  splitPos <- NULL
  
  # Split by ; to get the dp info
  splitPos <- strsplit(qualityVectorP[index], ";")[[1]][1]
  
  # Then split by = to get just the value
  splitPos <- strsplit(splitPos, "=")[[1]][2]
  
  # Convert to int
  splitPos <- as.integer(splitPos)
  
  # Insert the value into the dpVector
  dpVectorP[index] <- splitPos
  
  
}


# Get the mean
meanPDP <- mean(dpVectorP)

# Get the standard deviation
sdPDP <- sd(dpVectorP)

# Get the 95th percentile
the95thP <- qnorm(0.95, meanPDP, sdPDP)

# Convert the paired vectors back to a data frame
positionDPP <- data.frame(positionVectorP, dpVectorP)

# Get data for positions with depths above 140
above225P <- plotAtDepth(positionDPP, 225)

# Retrieve the regions
regions225P <- getMinMaxRegions(above225P, 500)

# Fix the gff table
gffTable <- processInfo(gffTable)

# Remove dbxref rows
gffTable = gffTable[-grep("Dbxref", gffTable$V9),]

# Get the start and end of full regions
regions225P <- getAnnoInfo(regions225P, gffTable)# Save plot as .pdf file

# Write to pdf
outputFile <- "CITP-Potential-Duplications.pdf"
pdf(outputFile, height=15, width=10)

# Create the plots
plotRegionCoverage(regions225P, positionDPP, gffTable)

# Close map
dev.off()

# Write things to files
write.xlsx(regions225P, "CITP-Potential-Duplicates.xlsx")
#############
##FUNCTIONS##
#############

# Function to plot positions above given read depth, and return df above depth
plotAtDepth <- function(dFrame, depth){
  
  # Create dataframe with reads only above depth
  aboveDepth <- subset.data.frame(dFrame, dFrame$dpVector >= depth)
  
  # Plot the positions
  plot(x = aboveDepth$positionVector, y = aboveDepth$dpVector,
       xlab = "Position", ylab = "Read Depth", type = "p",
       main = paste("Positions w/ Depth >=", toString(depth), sep = " "),
       pch=19, col=rgb(0,0,0, 0.1))
  
  plot(x = aboveDepth$positionVector, y=1:nrow(aboveDepth),
       xlab = "Position", ylab = "Index", type = "p",
       main = paste("Positions w/ Depth >=", toString(depth), sep = " "),
       pch=19, col=rgb(0,0,0, 0.1))
  
  return(aboveDepth)
  
}

# Function to retrieve regions from input following a distance threshold
getMinMaxRegions <- function(dFrame, threshold){
  
  # Set the first value of dFrame to equal the min, initiate vectors for mins/maxs
  currmin = dFrame[1, 1]
  currmax = dFrame[2, 1]
  minvect = c()
  maxvect = c()
  
  # Loop thru the dFrame and check where region ends
  for(row in 3:nrow(dFrame)){
    
    # Check to see if the current max is within threshold of the previous
    if(dFrame[row, 1] > (currmax + threshold)){
      
      # If outside, store the mins and maxes in vectors
      minvect = append(minvect, currmin)
      maxvect = append(maxvect, currmax)
      
      # Set value as the new min and max
      currmin = dFrame[row, 1]
      currmax = dFrame[row, 1]
      
    } else {
      
      # Set value as new max
      currmax = dFrame[row, 1]
    }
    
  }
  
  # Don't forget to write the last min/max to vector
  minvect = append(minvect, currmin)
  maxvect = append(maxvect, currmax)
  
  # Now place these into a dataframe
  minmaxFrame = data.frame(minvect, maxvect)
  
  # Create empty columns
  minmaxFrame[,c(3,4,5)] <- NA
  
  return(minmaxFrame)
}

# Function to process informations in gff table
processInfo <- function(annotable){
  
  # Loop thru each row 
  for(row in 1:nrow(annotable)){
    
    # split the info column
    splitter = strsplit(annotable[row, 3], split = ";")[[1]]
    
    # Store chosen info back in table
    annotable[row, 3] = paste(splitter[2], splitter[7], sep = ";")
  }
  
  return(annotable)
}

# Function to retrieve annotation info for each start/end
getAnnoInfo <- function(startend, annotable){
  
  # Loop thru each row in the startend frame
  for(row in 1:nrow(startend)){
    
    # Find where in the annotations the positions are
    location1 = length(which(annotable$V4 <= startend[row, 1]))
    location2 = which(annotable$V5 >= startend[row, 2])[1]
    
    # If they match, the location is in 1 gene, add to matrix
    if(location1 == location2){
      
      startend[row,3] = annotable[location1,1]
      startend[row,4] = annotable[location1,2]
      startend[row,5] = annotable[location1,3]
    
    # Else it's 2 genes, so store both of them
    } else{
      
      # Store bits and paste
      startend[row, 3] = annotable[location1,1]
      startend[row, 4] = annotable[location2,2]
      startend[row, 5] = paste(annotable[location1,3],annotable[location2,3],
                               sep = ";")
    }
    
  }
  
  return(startend)
}

# Function to plot coverage graphs for regions
plotRegionCoverage <- function(startend, dFrame, annotable){
  
  # Loop thru each row in the start/end frame
  for(row in 1:nrow(startend)){
    
    # Create variables to store from and to sites
    from = startend[row, 3]
    to = startend[row, 4]
    
    # Use the start/end position to access indexes of full frame and plot
    plot(x = dFrame[from:to, 1], y = dFrame[from:to, 2], 
         main = paste("Read depth for positions", toString(from), "to", toString(to)),
         xlab = "Position", ylab = "Read Depth", type = "l")
    
    # Use the from and tos to get the genes that they correspond to in the gff
    location = which(annotable$V4 == from)
    
    # Check if the plot is just this gene
    if(to == annotable[location, 2]){
      
      # Plot it 
      lines(x = dFrame[from:to, 1], y = dFrame[from:to, 2], col = "red")
      
      # Split off the relevant info
      splitter = strsplit(annotable[location, 3], split = ";")[[1]][1]
      splitter = strsplit(splitter, split = "-")[[1]][2]
      
      # Add legend
      legend("topright", 
             legend = splitter,
             text.col = "red", cex = 0.6)
      
    } else{
      
      # Assign actual tos and froms according to genes
      to1 = annotable[location, 2]
      from2 = annotable[location+1, 1]
      
      # Plot first gene
      lines(x = dFrame[from:to1, 1], y = dFrame[from:to1, 2], col = "red")
      
      # Plot second
      lines(x = dFrame[from2:to, 1], y = dFrame[from2:to, 2], col = "blue")
      
      # Split off the relevant info
      splitter = strsplit(annotable[location, 3], split = ";")[[1]][1]
      splitter = strsplit(splitter, split = "-")[[1]][2]
      
      splitter2 = strsplit(annotable[location+1, 3], split = ";")[[1]][1]
      splitter2 = strsplit(splitter2, split = "-")[[1]][2]
      
      # Add legend
      legend("topright", 
             legend = c(splitter, splitter2, "Intergenic region"),
             text.col = c("red", "blue", "black"), cex = 0.6)
    }
    
  }
}



