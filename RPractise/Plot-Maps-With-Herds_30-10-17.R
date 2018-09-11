# Set path variable
pathCoords <- "path/to/polygons" #general path for all files
pathSamples <- "path/to/csv" #samples path

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Londonderry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Read in table
fileName <- paste(pathSamples, "GenotypingData.csv", sep="")
sampleTable <- read.table(fileName,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                          check.names=FALSE) # Names left as they are, no dots inserted

# Initialise vector with column headers
colHeads <- c("292 Code", "X3 Code", "25 Code", "47 Code", "3 Code", "7 Code", 
              "10 Code",  "32 Code")

# Save plot as .pdf file
outputFile <- paste(pathSamples, "VNTRMaps.pdf", sep="")
pdf(outputFile, height=35, width=14)

# Configure plot layout
layout(matrix(c(1,1,1,1,2,3,4,5,6,7), nrow=5, ncol=2, byrow=TRUE))

# Count the number of samples per county
sampleNumbers <- numberSamples(sampleTable)

# Count number of samples per herd
sampleHerds <- herdsSamples(sampleTable)

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

# Plot county polygons and related sample data
polygonsData(polygonCoords, sampleNumbers, counties)

# Find unique VNTRs
VNTRUniques <- findUniqueVNTRS(sampleTable, colHeads)
UniqueNames <- c("INMV 2", "INMV 1", "NA", "No Match, but probably 1", "INMV 116",
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

numberSamples <- function(sampleTable) { 
  
  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()
  
  # Split off the county name away from the other contents in the herd location strings
  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available
    if(is.na(sampleTable[row, "292 Code"]) == TRUE){
      next
    }
    
    # Get county name from current row
    # Split away the second part of a string for every row in the herd loc. column
    countyName <- strsplit(sampleTable[row, "Herd Location"], split = " ")[[1]][1] 
    
    # Check if we have encountered the current county before
    if(is.null(countySamples[[countyName]]) == TRUE) { # if county name not there
      countySamples[[countyName]] <- 1 # add it into the list and count once
    } else { # if already there add an extra 1 to the count for each time seen
      countySamples[[countyName]] <- countySamples[[countyName]] + 1 
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
    if(is.na(sampleTable[row, "292 Code"]) == TRUE){
      next
    }
    
    # Get county name from current row
    # Split away the second part of a string for every row in the herd loc. column
    countyName <- strsplit(sampleTable[row, "Herd Location"], split = " ")[[1]][1] 
    
    # Get VNTR type
    pasteVNTR <- paste(sampleTable[row, colHeads], collapse = ":") 
    
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
    
    # Convert VNTRs for each county into a data frame to get frequency
    VNTRPresence[[index]] <- as.data.frame(table(VNTRPresence[[keys[index]]]))
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
      }else if(type %in% as.character(frequentVNTRs[[county]][,1]) == FALSE){
        
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
        typeRow <- which(as.character(frequentVNTRs[[county]][,1]) == type)
        
        # Get number of samples for current type
        nSamples <- frequentVNTRs[[county]][typeRow,2]
        
        # Calculate proportion of total number of samples for current type
        proportion <- nSamples/VNTROverall[[type]]
        
        # Plot the polygon with county name, number of samples in brackets
        # Colour polygon by proportion
        polygon(x=polygonCoords[[county]][, "X"],
                y=polygonCoords[[county]][, "Y"],
                col=rgb(red=0, green=0, blue=1, alpha=proportion)) 
        
        # Build label by pasting strings together
        countyLabel <- paste(county, " ",
                             "(",nSamples, ")", sep = "")
        
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
    if(sampleNumbers[[keys[index]]] > currentmaxValue) { # if they are greater than 1
      currentmaxValue = sampleNumbers[[keys[index]]] # current max is changed to that
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
      nSamples <- sampleNumbers[[counties[index]]]
      
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
                           "(",nSamples, ")", sep = "")
      
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

### Herd Functions ###

herdsSamples <- function(sampleTable) { 
  
  # Initialise an empty list to store herd names and sample count for each herd
  herdSamples <- list()
  

  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available
    if(is.na(sampleTable[row, "292 Code"]) == TRUE){
      next
    }
    
    # Get herd name from current row
    herdName <- sampleTable[row, "Herd Location"] 
    
    # Check if we have encountered the current herd before
    if(is.null(herdSamples[[herdName]]) == TRUE) { # if herd name not there
      herdSamples[[herdName]] <- 1 # add it into the list and count once
    } else { # if already there add an extra 1 to the count for each time seen
      herdSamples[[herdName]] <- herdSamples[[herdName]] + 1 
    }
  }
  return(herdSamples)
}

herdsVNTRS <- function(sampleTable, colHeads) { 
  
  # Initialise an empty list to store county names and VNTRs for each county
  herdVNTRs <- list()
  
  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available
    if(is.na(sampleTable[row, "292 Code"]) == TRUE){
      next
    }
    
    # Get herd name from current row
    herdName <- sampleTable[row, "Herd Location"] 
    
    # Get VNTR type
    pasteVNTR <- paste(sampleTable[row, colHeads], collapse = ":") 
    
    emptyVector <- c()
    
    # Check whether you have already added current herd in
    if(is.null(herdVNTRs[[herdName]]) == FALSE){
      
      # Get VNTR types of current herd
      emptyVector <- herdVNTRs[[herdName]]
      # Append current VNTR into vector
      emptyVector[length(emptyVector) + 1] <- pasteVNTR
      # Store vector in list under current county
      herdVNTRs[[herdName]] <- emptyVector
    }else{
      
      # Store current VNTR in vector in list under herd name as key
      emptyVector <- pasteVNTR
      herdVNTRs[[herdName]] <- emptyVector  
    }
    
  }
  return(herdVNTRs)
}

freqherdVNTR <- function(VNTRHerdPresence) {
  
  # Set herd names as keys to list
  keys <- names(VNTRHerdPresence)
  
  # For each element of the list
  for(index in 1:length(keys)) {
    
    # Convert VNTRs for each herd into a data frame to get frequency
    VNTRHerdPresence[[index]] <- as.data.frame(table(VNTRHerdPresence[[keys[index]]]))
  }
  return(VNTRHerdPresence)
  
}