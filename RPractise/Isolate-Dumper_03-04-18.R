# Set path variables
pathIso <- "path/to/csv"
pathCov <- "path/to/coverage/summary"
pathExtra <- "path/to/list"

# Read in table of isolates
isoTable <- read.table(pathIso,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                          check.names=FALSE) # Names left as they are, no dots inserted

# Read in table of isolate coverage
covTable <- read.table(pathCov,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted

# Read in table of additional isolates to dump
extraTable <- read.table(pathExtra,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                         check.names=FALSE) # Names left as they are, no dots inserted



# Run function to generate file
checkDup(isoTable, covTable, extraTable)

# Function to check for duplication of isolates and pick those with best coverage
checkDup <- function(isoTable, covTable, extraTable){
  
  # Create vector to store dump isolates
  dumpIsolates <- c()
  
  # Loop thru rows of isolate table
  for(row in 1:nrow(isoTable)){
    
    # If nothing there in secondary accession, skip row
    if(isoTable[row, "Secondary Accession"] == ""){
      next
    } else{
      
      # Check where in the coverage vector the Primary and Secondary accessions are
      locationPri <- grep(isoTable[row, "Accession"], covTable[,1])
      locationSec <- grep(isoTable[row, "Secondary Accession"], covTable[,1])
      
      # Compare the coverage values of the two, place the worst one in dump list
      if(covTable[locationPri, 3] > covTable[locationSec, 3]){
        
        dumpIsolates <- append(dumpIsolates, covTable[locationSec, 1])
      } else {
        
        dumpIsolates <- append(dumpIsolates, covTable[locationPri, 1])
      } 
    }
  }
  
  # Loop thru each row in the extra dump table 
  for(row in 1:nrow(extraTable)){
    
    # Check where in the coverage vector the accessions are
    locationAcc <- grep(extraTable[row, "Access"], covTable[,1])
    
    # Dump the names into the dump vector
    dumpIsolates <- append(dumpIsolates, covTable[locationAcc, 1])
    
  }
  
  theDump <- data.frame(dumpIsolates)
    
  # Write the dumped isolates to file
  cat(capture.output(theDump), file = "isolate_to_dump.txt", sep = "\n")
  
}