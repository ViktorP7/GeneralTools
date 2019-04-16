
# Set the path
path <- "path/to/csv"

numberSamples(path)

numberSamples <- function(path) { 

  # Initialise an empty list to store county names and sample count for each county
  countySamples <- list()

  # Read in data from .csv file
  fileName <- paste(path, "GenotypingData.csv", sep="")
  sampleTable <- read.table(fileName,
             header = TRUE,
             sep = ",",
             stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
             check.names=FALSE) # Names left as they are, no dots inserted

  # Split off the county name away from the other contents in the herd location strings
  for(row in 1:nrow(sampleTable)) { # for each row
  
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
