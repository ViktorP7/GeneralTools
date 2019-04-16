# Set path variable
pathCoords <- "path/to/polygons" #general path for all files
pathSamples <- "path/to/samples" #samples path

# Count the number of samples per county
sampleNumbers <- numberSamples(pathSamples, "GenotypingData.csv")

# Plot the polygons of the following counties on the same figure
# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
  "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
  "Limerick","Londonderry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
  "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

polygonCoords <- getPolygonCoords(counties, pathCoords)

# Calculate limits of plot
ranges <- mapLimits(polygonCoords, counties)

# Calculate the max number of samples
maxNSamples <- sampleMax(sampleNumbers)

# Save plot as .pdf file
outputFile <- paste(pathSamples, "sampleMap.pdf", sep="")
pdf(outputFile, height=7, width=7)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(ranges[1], ranges[2]), 
     ylim = c(ranges[3], ranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Plot county polygons and related sample data
polygonsData(polygonCoords, sampleNumbers, counties)

# Save map
dev.off()

# Check if key is present in list
#is.null(list[[key]])

# Rule for functions
# Everything referred to in function must be passed in or created in function

#############
# FUNCTIONS #
#############

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
  
addPolygons <- function(counties, polygonCoords) {
    
  # Add each polygon to empty plot
  for(index in 1:length(counties)) {
    polygon(x=polygonCoords[[counties[index]]][, "X"],
            y=polygonCoords[[counties[index]]][, "Y"])
    
    # Find rough centre of county polygon and add appropriate names to each county centre
    text(x=mean(polygonCoords[[counties[index]]][,"X"]),
         y=mean(polygonCoords[[counties[index]]][,"Y"]),
         labels = counties[index])
  }
}  

numberSamples <- function(path, fileName) { 
  
  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()
  
  # Read in data from .csv file
  fileName <- paste(path, fileName, sep="")
  sampleTable <- read.table(fileName,
                            header = TRUE,
                            sep = ",",
                            stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                            check.names=FALSE) # Names left as they are, no dots inserted
  
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
           cex = 0.7) # sets font size
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
           cex = 0.7) 
    }
  }
}