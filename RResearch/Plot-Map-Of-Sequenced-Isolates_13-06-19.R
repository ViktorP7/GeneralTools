# 13-06-19 Script created to work side by side with "Plot True Irish Tree" script

# Set path variables
pathCoords <- "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Find the distances between all isolates
allDist <- cophenetic(irishOnlytree)

# Pull out herd info from tip labels
allHerds <- isolateHerder(irishOnlytree$tip.label)

# Get herd distances
herdDist <- allDist[allHerds, allHerds]

# Process distance matrices for each herd
herdDist[upper.tri(herdDist)] <- NA

# Get the herd names
herdNames <- getHerdNames(herdDist)

# Count the number of samples per county
sampleNumbers <- numberSamples(herdNames)

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
legend("topleft", legend="(N. isolates/N. herds)", bty="n", cex=1.8)

# Plot county polygons and related sample data
polygonsData(polygonCoords, sampleNumbers, counties)

# Plot small map
smallMap(polygonCoords, counties)

#############
# Functions #
#############

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
    
    # Check what colour needs to be assigned 
    if(counties[index] == "Donegal" ||counties[index] == "Derry" 
       ||counties[index] == "Monaghan" ||counties[index] == "Cavan"||counties[index] == "Fermanagh"
       ||counties[index] == "Tyrone"||counties[index] == "Armagh"||counties[index] == "Down"
       ||counties[index] == "Antrim"){
      
      colour = "black"
    } else if(counties[index] == "Mayo" || counties[index] =="Roscommon" 
              ||counties[index] == "Sligo" ||counties[index] == "Leitrim"
              ||counties[index] == "Galway"){
      
      colour = "deepskyblue4"
    } else if(counties[index] == "Clare" ||counties[index] == "Kerry" 
              ||counties[index] == "Cork" ||counties[index] == "Limerick" 
              ||counties[index] == "Tipperary" ||counties[index] == "Waterford"){
      
      colour = "darkorange4"
    } else {
      
      colour = "red"
    }
    
    # Find rough centre of county polygon and add appropriate names to each county centre
    meanX <- mean(polygonCoords[[counties[index]]][,"X"])
    meanY <- mean(polygonCoords[[counties[index]]][,"Y"])
    
    # Add county name
    # colour by whether or not sampled - red = sampled, grey = unsampled
    # If the county name is not present in the sample list, it's grey
    if(is.null(sampleNumbers[[counties[index]]]) == TRUE ){
      
      # Plot the polygon of the current county
      polygon(x=polygonCoords[[counties[index]]][, "X"],
              y=polygonCoords[[counties[index]]][, "Y"],
              border = colour)
      
      # Add county name
      text(x= meanX,
           y= meanY,
           labels = counties[index],
           col = "grey",
           cex = 1.0) # sets font size
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
              col=rgb(red=0, green=0, blue=1, alpha=proportion),
              border = colour) 
      
      # Build label by pasting strings together
      countyLabel <- paste(counties[index], " ",
                           "(",nSamples, "/", nHerdsSamples,
                           ")", sep = "")
      
      # Add label
      text(x= meanX,
           y= meanY,
           labels = countyLabel,
           col = colour,
           cex = 1.0) 
    }
  }
}

smallMap <- function(polygonCoords, counties) {
  
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
