# Set path variables
path <- "path/to/directory"
pathCoords <- paste(path, "PolygonCoords/", sep="") #general path for all files

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Londonderry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Read in table
fileName <- paste(path, "AllGenosmod.csv", sep="")
sampleTable <- read.table(fileName,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                          check.names=FALSE) # Names left as they are, no dots inserted

# Get the isolates
ourIsolates <- getIsolates(path, "corkisolates.txt")

# Initialise vector with column headers
colHeads <- c("292 Code", "X3 Code", "25 Code", "47 Code", "3 Code", "7 Code", 
              "10 Code",  "32 Code")

# Save plot as .pdf file
outputFile <- paste(path, "AllMAPSampling_14-06-18.pdf", sep="")
pdf(outputFile, height=35, width=14)

# Configure plot layout
layout(matrix(c(1,1,1,1,2,3,4,5,6,7), nrow=5, ncol=2, byrow=TRUE))

# Edit the margin sizes
par(mar=c(1,1,1,1)) # Bottom, left, top, right

# Count the number of samples per county
sampleNumbers <- numberSamples(sampleTable, ourIsolates)

# Get polygon coordinates for map plotting
polygonCoords <- getPolygonCoords(counties, pathCoords)

# Calculate limits of plot
ranges <- mapLimits(polygonCoords, counties)

# Calculate the max number of samples
maxNSamples <- sampleMax(sampleNumbers)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Add legend
legend("topleft", legend="(N. isolates/N. herds)", bty="n", cex=3)

# Plot county polygons and related sample data
polygonsData(polygonCoords, sampleNumbers, counties)

# Find unique VNTRs
VNTRUniques <- findUniqueVNTRS(sampleTable, colHeads)
UniqueNames <- c("INMV 2", "INMV 1", "NA", "No Match, could be 1", "INMV 116",
                 "No Match, could be 1", "INMV 3")
# Place Unique VNTRs in a list
VNTRIdentifiers <- UniquesList(VNTRUniques, UniqueNames) 

# Find VNTRs present in each county
VNTRPresence <- numberVNTRS(sampleTable,colHeads)

# Get freqeuncy of VNTRs in counties
frequentVNTRs <- freqVNTR(VNTRPresence)

# Count overall presence of VNTRs
VNTROverall <- overallVNTRTypes(sampleTable, colHeads)

# Make maps
mapping(VNTRUniques, pathSamples, counties, polygonCoords, frequentVNTRs, 
        ranges, VNTRIdentifiers)

# Close map
dev.off()

#par( mfrow = c( 2, 3 ) )


#pdf(fileName, height=35, width=14)

#table <- table[, c(2,3,1,5,4)]

#############
# Functions #
#############

numberSamples <- function(sampleTable, ourIsolates) { 
  
  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()
  countyHerds <- list()
  
  # Split off the county name away from the other contents in the herd location strings
  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available, check if isolate is present
    if(is.na(sampleTable[row, "292 Code"]) == TRUE || 
       sampleTable[row, "MAP Isolates (TB14-00)"] %in% ourIsolates == FALSE){
      next
    }
    
    # Get county name from current row
    # Split away the second part of a string for every row in the herd loc. column
    nameParts <- strsplit(sampleTable[row, "Herd Location"], split = " ")[[1]]
    herdSampled <- nameParts[2]
    countyName <- nameParts[1]
    
    # Check if we have encountered the current county before
    if(is.null(countySamples[[countyName]]) == TRUE) { # if county name not there
      countySamples[[countyName]] <- c(1, 1) # add it into the list and count once
      countyHerds[[countyName]] <- c(herdSampled)
    } else { # if already there add an extra 1 to the count for each time seen
      
      countySamples[[countyName]][1] <- countySamples[[countyName]][1] + 1
      
      if(herdSampled %in% countyHerds[[countyName]] == FALSE){
        countyHerds[[countyName]] <- c(countyHerds[[countyName]], herdSampled)
        countySamples[[countyName]][2] <- length(countyHerds[[countyName]])
      }
    }
  }
  return(countySamples)
}

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

findUniqueVNTRS <- function(sampleTable, colHeads) {
  
  # Initialise vector to store VNTR types
  vntrTypes <- c()
  
  # For each row in all rows in table
  for(row in 1:nrow(sampleTable)) {
    vntrTypes[row] <- paste(sampleTable[row, colHeads], collapse = ":")
  }
  
  vNTRs <- unique(vntrTypes)
  
  #vntrTypes <- vntrTypes[grepl(vntrTypes, pattern="NA") == FALSE]
  #barplot(table(vntrTypes), horiz=TRUE, las=1)
  
  return(vNTRs)
}

numberVNTRS <- function(sampleTable, colHeads) { 
  
  # Initialise an empty list to store county names and VNTRs for each county
  countyVNTRs <- list()
  
  # Split off the county name away from the other contents in the herd location strings
  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available
    # Check if isolate is present
    if(is.na(sampleTable[row, "292 Code"]) == TRUE || 
       sampleTable[row, "MAP Isolates (TB14-00)"] %in% ourIsolates == FALSE){
      next
    }
    
    # Get county name from current row
    # Split away the second part of a string for every row in the herd loc. column
    nameParts <- strsplit(sampleTable[row, "Herd Location"], split = " ")[[1]]
    countyName <- nameParts[1]
    herdSampled <- nameParts[2]
    
    # Get VNTR type
    pasteVNTR <- paste(sampleTable[row, colHeads], collapse = ":") 
    
    # Add on Herd
    pasteVNTR <- paste(pasteVNTR, herdSampled, sep="-")
    
    emptyVector <- c()
    
    # Check whether you have already added current county in
    if(is.null(countyVNTRs[[countyName]]) == FALSE){
      
      # Get VNTR types of current county
      emptyVector <- countyVNTRs[[countyName]]
      # Append current VNTR into vector
      emptyVector[length(emptyVector) + 1] <- pasteVNTR
      # Store vector in list under current county
      countyVNTRs[[countyName]] <- emptyVector
    }else{
      
      # Store current VNTR in vector in list under county name as key
      emptyVector <- pasteVNTR
      countyVNTRs[[countyName]] <- emptyVector  
    }
    
  }
  return(countyVNTRs)
}

freqVNTR <- function(VNTRPresence) {
  
  # Set county names as keys to list
  keys <- names(VNTRPresence)
  
  # For each element of the list
  for(index in 1:length(keys)) {
    
    # Get the VNTR-Herd information
    vntrHerdInfoForCounty <- VNTRPresence[[keys[index]]]
    
    # Initialise table with N. samples and N. herds sampled per VNTR type
    vntrInfo <- data.frame(VNTR=c(NA), nSamples=c(NA), nHerdsSampled=c(NA),
                           stringsAsFactors=FALSE)
    
    # Initialise a list that records herds for VNTRs
    vntrHerds <- list()
    
    # Examine VNTR-herd info
    for(i in 1:length(vntrHerdInfoForCounty)){
      
      # Get the VNTR
      parts <- strsplit(vntrHerdInfoForCounty[i], split="-")[[1]]
      vntr <- parts[1]
      herd <- parts[2]
      
      # Check if encountered current VNTR
      if(vntr %in% vntrInfo[, "VNTR"] == TRUE){
        
        # Find row of VNTR
        row <- which(vntrInfo[, "VNTR"] == vntr)
        
        # Update count for current VNTR
        vntrInfo[row, "nSamples"] <- vntrInfo[row, "nSamples"] + 1
        
        # Check if encountered herd for this vntr
        if(herd %in% vntrHerds[[vntr]] == FALSE){
          vntrHerds[[vntr]] <- c(vntrHerds[[vntr]], herd)
          
          # Update number of herds associated with vntr
          vntrInfo[row, "nHerdsSampled"] <- length(vntrHerds[[vntr]])
        }
        
      }else{
        vntrInfo[nrow(vntrInfo) + 1, "VNTR"] <- vntr
        vntrInfo[nrow(vntrInfo), "nSamples"] <- 1
        vntrInfo[nrow(vntrInfo), "nHerdsSampled"] <- 1
      }
    }
    
    # Convert VNTRs for each county into a data frame to get frequency
    VNTRPresence[[index]] <- vntrInfo[-1, ]
  }
  return(VNTRPresence)
  
}

mapping <- function(VNTRUniques, pathSamples, counties, polygonCoords, frequentVNTRs,
                    ranges, VNTRIdentifiers){
  
  # Plot each VNTR onto a separate map
  for(type in VNTRUniques){
    if(type == "NA:NA:NA:NA:NA:NA:NA:NA"){
      next
    }
    
    # Create empty plot with correct dimensions
    plot(x=NA, y=NA,
         xlim = c(ranges[1], ranges[2]), 
         ylim = c(ranges[3], ranges[4]),
         main = paste(type, VNTRIdentifiers[[type]]), xlab = "", ylab = "",
         bty = "n", axes = FALSE)
    
    # Examine sample counts for each county
    for(county in counties) {
      
      # Find rough centre of county polygon and add appropriate names to each county centre
      meanX <- mean(polygonCoords[[county]][,"X"])
      meanY <- mean(polygonCoords[[county]][,"Y"])
      
      # Plot polygons of completely unsampled counties
      if(is.null(frequentVNTRs[[county]]) == TRUE ){
        
        # Plot the polygon of the current county
        polygon(x=polygonCoords[[county]][, "X"],
                y=polygonCoords[[county]][, "Y"])
        
        # Add county name
        text(x= meanX,
             y= meanY,
             labels = county,
             col = "grey",
             cex = 0.7) # sets font size
        
        # Plot polygons of counties with no samples of current VNTR type
      }else if(type %in% as.character(frequentVNTRs[[county]][, "VNTR"]) == FALSE){
        
        # Plot the polygon of the current county
        polygon(x=polygonCoords[[county]][, "X"],
                y=polygonCoords[[county]][, "Y"])
        
        # Add county name
        text(x= meanX,
             y= meanY,
             labels = county,
             col = "purple",
             cex = 0.7) # sets font size
        
        # Plot polygons of counties with samples of current VNTR type 
      }else{
        
        # Find the row that the current type is in in the current counties frequency table
        typeRow <- which(frequentVNTRs[[county]][, "VNTR"] == type)
        
        # Get number of samples for current type
        nSamples <- frequentVNTRs[[county]][typeRow, "nSamples"]
        nHerdsSampled <- frequentVNTRs[[county]][typeRow, "nHerdsSampled"]
        
        # Calculate proportion of total number of samples for current type
        proportion <- nSamples/VNTROverall[[type]]
        
        # Plot the polygon with county name, number of samples in brackets
        # Colour polygon by proportion
        polygon(x=polygonCoords[[county]][, "X"],
                y=polygonCoords[[county]][, "Y"],
                col=rgb(red=0, green=0, blue=1, alpha=proportion)) 
        
        # Build label by pasting strings together
        countyLabel <- paste(county, " ",
                             "(", nSamples, "/", nHerdsSampled, ")", sep = "")
        
        # Add label
        text(x= meanX,
             y= meanY,
             labels = countyLabel,
             col = "red",
             cex = 0.7)
        
      }
    }
  }
}

overallVNTRTypes <- function(sampleTable, colHeads) { 
  
  # Initialise an empty list to store VNTR counts for each VNTR type
  vNTRCount <- list()
  
  # For each row in the table
  for(row in 1:nrow(sampleTable)) { 
    
    # Paste together VNTR columns to get VNTR type
    vNTRPaste <- paste(sampleTable[row, colHeads], collapse = ":") 
    
    # Check for NaN VNTR types
    if(vNTRPaste == "NA:NA:NA:NA:NA:NA:NA:NA"){
      next # Skips to next iteration of for loop, or can use break to exit loop
    }
    
    # Check if we have encountered the current VNTR before
    if(is.null(vNTRCount[[vNTRPaste]]) == TRUE) { # if VNTR not there
      vNTRCount[[vNTRPaste]] <- 1 # add it into the list and count once
    } else { # if already there add an extra 1 to the count for each time seen
      vNTRCount[[vNTRPaste]] <- vNTRCount[[vNTRPaste]] + 1 
    }
  }
  return(vNTRCount)
} 

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

polygonsData <- function(polygonCoords, sampleNumbers, counties) {
  
  for(index in 1:length(counties)) {
    
    # Find rough centre of county polygon and add appropriate names to each county centre
    meanX <- mean(polygonCoords[[counties[index]]][,"X"])
    meanY <- mean(polygonCoords[[counties[index]]][,"Y"])
    
    # Add county name
    # colour by whether or not sampled - red = sampled, grey = unsampled
    # If the county name is not present in the sample list, it's grey
    if(is.null(sampleNumbers[[counties[index]]]) == TRUE ){
      
      # Plot the polygon of the current county
      polygon(x=polygonCoords[[counties[index]]][, "X"],
              y=polygonCoords[[counties[index]]][, "Y"])
      
      # Add county name
      text(x= meanX,
           y= meanY,
           labels = counties[index],
           col = "grey",
           cex = 1.2) # sets font size
    }else{ # counties that are present on the sample list
      
      # Get number of samples for current county
      nSamples <- sampleNumbers[[counties[index]]][1]
      nHerdsSamples <- sampleNumbers[[counties[index]]][2]
      
      # Calculate alpha = proportion of max number samples at current county
      proportion <- nSamples/maxNSamples
      
      # Plot the polygon of the current county
      # RGB allows for any colour, alpha is the transparency of the colour
      # RGB values go from 0 to 1
      polygon(x=polygonCoords[[counties[index]]][, "X"],
              y=polygonCoords[[counties[index]]][, "Y"],
              col=rgb(red=0, green=0, blue=1, alpha=proportion)) 
      
      # Build label by pasting strings together
      countyLabel <- paste(counties[index], " ",
                           "(",nSamples, "/", nHerdsSamples,
                           ")", sep = "")
      
      # Add label
      text(x= meanX,
           y= meanY,
           labels = countyLabel,
           col = "red",
           cex = 1.2) 
    }
  }
}

UniquesList <- function(VNTRUniques, UniqueNames) {
  VNTRIDs <- list()
  
  for(index in 1:length(VNTRUniques)) {
    VNTRIDs[[VNTRUniques[index]]] <- UniqueNames[index]
  }
  return(VNTRIDs) 
} 

# Function to get isolate codes
getIsolates <- function(path, nameOfFile){
  
  # Read in the sample texts
  fileNameTwo <- paste(path, nameOfFile, sep = "")
  sampleLines <- readLines(fileNameTwo, warn = FALSE)
  
  # Split off the TB14/13 and (P) from the isolates
  # New vector to store isolates
  isolateNos <- c()
  
  # Loop thru isolates and split each one to just the number
  for(isolate in sampleLines){
    
    curriso <- ""
    
    curriso <- strsplit(isolate, " ")[[1]][2]
    
    isolateNos <- append(isolateNos, curriso)
  }
  
  return(isolateNos)
}
