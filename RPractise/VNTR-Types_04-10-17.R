# Define path to table
pathSamples <- "path/to/files"

# Read in table
fileName <- paste(pathSamples, "GenotypingData.csv", sep="")
sampleTable <- read.table(fileName,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                          check.names=FALSE) # Names left as they are, no dots inserted

# Initialise vector with column headers
colHeads <- c("X3 Code", "3 Code", "7 Code", "10 Code", "25 Code", "32 Code",
              "47 Code", "292 Code")

# Find unique VNTRs
VNTRUniques <- findUniqueVNTRS(sampleTable, colHeads)

# Find VNTRs present in each county
VNTRPresence <- numberVNTRS(sampleTable,colHeads)

# Count overall presence of VNTRs
VNTROverall <- overallVNTRTypes(sampleTable, colHeads)

# Get freqeuncy of VNTRs in counties
frequentVNTRs <- freqVNTR(VNTRPresence)

#############
# Functions #
#############
findUniqueVNTRS <- function(sampleTable, colHeads) {
  
  # Initialise vector to store VNTR types
  vntrTypes <- c()
  
  # For each row in all rows in table
  for(row in 1:nrow(sampleTable)) {
    vntrTypes[row] <- paste(sampleTable[row, colHeads], collapse = ":")
  }

  vNTRs <- unique(vntrTypes)
  
  vntrTypes <- vntrTypes[grepl(vntrTypes, pattern="NA") == FALSE]
  barplot(table(vntrTypes), horiz=TRUE, las=1)

  return(vNTRs)
}  

numberVNTRS <- function(sampleTable, colHeads) { 
  
  # Initialise an empty list to store county names and VNTRs for each county
  countyVNTRs <- list()
  
  # Split off the county name away from the other contents in the herd location strings
  for(row in 1:nrow(sampleTable)) { # for each row
    
    # Check if VNTR information available
    if(is.na(sampleTable[row, "X3 Code"]) == TRUE){
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

overallVNTRTypes <- function(sampleTable, colHeads) { 
    
    # Initialise an empty list to store VNTR counts for each VNTR type
    vNTRCount <- list()
    
    # For each row in the table
    for(row in 1:nrow(sampleTable)) { 
      
      # Paste together VNTR columns to get VNTR type
      vNTRPaste <- paste(sampleTable[row, colHeads], collapse = ":") 
      
      # Check for NaN VNTR types
      if(vNTRPaste == "NANANANANANANANA"){
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
