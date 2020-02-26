# 13-06-19 Script created to work side by side with "Plot True Irish Tree" script
# 12-08-19 Replaced obsolete functions with getNames function
# 19-08-19 Added arrows to map to show direction of movement between cattle birth/current counties
# 04-09-19 Arrows now weighted based on how many moves occur in the same way

# Load required packages
library(scales)

# Set path variables
pathCoords <- "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Get herd names from allDist
herdNames <- getNames(allDist, "Herd")

# Count the number of samples per county
sampleNumbers <- numberSamples(herdNames)

# Get polygon coordinates for map plotting
normalPolygonCoords <- getPolygonCoords(counties, pathCoords)

# Calculate limits of plot
normalranges <- mapLimits(normalPolygonCoords, counties)

# Calculate the max number of samples
maxNSamples <- sampleMax(sampleNumbers)

# Get the birth and current county of isolates
countyNames <- getNames(allDist, "CCounty")
birthcountyNames <- getNames(allDist, "BCounty")

# Save plot as .pdf file
outputFile <- paste("IrishArrowMap.pdf", sep="")
pdf(outputFile, height=75, width=75)

par(bg=NA)

# Create empty plot, with input of limits from above function
plot(x=NA, y=NA,
     xlim = c(normalranges[1], normalranges[2]), 
     ylim = c(normalranges[3], normalranges[4]),
     main = "", xlab = "", ylab = "",
     bty = "n", axes = FALSE)

# Plot county polygons and related sample data
polygonsData(normalPolygonCoords, sampleNumbers, counties)

# Add birth to current county arrows onto map
addArrows(countyNames, birthcountyNames, normalPolygonCoords)

dev.off()

# Plot small map
smallMap(polygonCoords, counties)

#############
# Functions #
#############

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
              border = alpha(colour, 0.6))
      
      # Add county name
      #text(x= meanX,
          # y= meanY,
           #labels = counties[index],
           #col = "grey",
           #cex = 1.0) # sets font size
    }else{ # counties that are present on the sample list
      
      # Get number of samples for current county
      nSamples <- sampleNumbers[[counties[index]]][1]
      nHerdsSamples <- sampleNumbers[[counties[index]]][2]
      
      # Calculate alpha = proportion of max number samples at current county
      #proportion <- nSamples/maxNSamples
      
      # Plot the polygon of the current county
      # RGB allows for any colour, alpha is the transparency of the colour
      # RGB values go from 0 to 1
      polygon(x=polygonCoords[[counties[index]]][, "X"],
              y=polygonCoords[[counties[index]]][, "Y"],
              #col=rgb(red=0, green=0, blue=1, alpha=proportion),
              border = colour) 
      
      # Build label by pasting strings together
      countyLabel <- paste(nSamples, "/", nHerdsSamples, sep = "")
      
      # Add label
      #text(x= meanX,
      #     y= meanY,
      #     labels = countyLabel,
      #     col = colour,
      #     cex = 2) 
    }
  }
}

addArrows <- function(countynames, birthcountynames, polygonCoords) {
  
  # Initialise vector to store birth/current combos
  combo <- c()
  
  # Loop thru each entry in the names
  for(index in 1:length(countynames)){
    
    # Check if there is an NA in the current position in either birth or current
    if(is.na(countynames[index]) || is.na(birthcountynames[index]) || birthcountynames[index] == "n/a"){
      
      next
    } else if(countynames[index] == birthcountynames[index]){
      
      next
    } else{
      
      # Append the combo into the combo vector
      combo <- append(combo, paste(birthcountynames[index], "_", countynames[index]))
    }
  }  
  
  # Create a vector of length combo
  dusted <- rep(NA, length(combo))
  
  # Loop thru the combo vector
  for(index in 1:length(combo)){
    
    if(combo[index] %in% dusted == TRUE){
      
      next
    } else{
      
      dusted[index] <- combo[index]
    
      # Weight arrow thickness based on how many occurences of a thing
      weight <- sum(combo == combo[index])
      
      # Split apart the names
      birth <- strsplit(combo[index], split = " _ ")[[1]][1]
      current <- strsplit(combo[index], split = " _ ")[[1]][2]
      
      # Get the rough centre of each polygon
      
      meanXB <- mean(polygonCoords[[birth]][,"X"])
      meanYB <- mean(polygonCoords[[birth]][,"Y"])
      meanXC <- mean(polygonCoords[[current]][,"X"])
      meanYC <- mean(polygonCoords[[current]][,"Y"])
      
      # Check what colour needs to be assigned 
      if(birth == "Donegal" ||birth == "Derry" 
          ||birth == "Monaghan" ||birth == "Cavan"
          ||birth == "Fermanagh" ||birth == "Tyrone"
          ||birth == "Armagh"||birth == "Down"
          ||birth == "Antrim"){
        
        colour = "black"
      } else if(birth == "Mayo" || birth =="Roscommon" 
                ||birth == "Sligo" ||birth == "Leitrim"
                ||birth == "Galway"){
        
        colour = "deepskyblue4"
      } else if(birth == "Clare" ||birth == "Kerry" 
                ||birth == "Cork" ||birth == "Limerick" 
                ||birth == "Tipperary" ||birth == "Waterford"){
        
        colour = "darkorange4"
      } else {
        
        colour = "red"
      }
      
      if(weight == 2){
        
        arrows(meanXB, meanYB, meanXC, meanYC, col = alpha(colour, 0.8), length = 0.1, lwd = weight+10)
      } else if(weight == 3){
        
        arrows(meanXB, meanYB, meanXC, meanYC, col = alpha(colour, 0.6), length = 0.1, lwd = weight+10)
      } else if(weight > 3){
        
        arrows(meanXB, meanYB, meanXC, meanYC, col = alpha(colour, 0.4), length = 0.1, lwd = weight+10)
      } else {
        
        arrows(meanXB, meanYB, meanXC, meanYC, col = colour, length = 0.1, lwd = weight+10)
      }
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
