
# Set the path
path <- "path/to/samples"

herdsSamples(path)

herdsSamples <- function(path) { 
  
  # Initialise an empty list to store herd names and sample count for each herd
  herdSamples <- list()
  
  # Read in data from .csv file
  fileName <- paste(path, "GenotypingData.csv", sep="")
  sampleTable <- read.table(fileName,
                            header = TRUE,
                            sep = ",",
                            stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                            check.names=FALSE) # Names left as they are, no dots inserted
  
  for(row in 1:nrow(sampleTable)) { # for each row
    
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
