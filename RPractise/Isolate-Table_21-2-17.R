# Set path variables
path <- "path/to/files"
pathCoords <- paste(path, "PolygonCoords/", sep="") #general path for coord files

# Initialise vector with column headers
colHeads <- c("292 Code", "X3 Code", "25 Code", "47 Code", "3 Code", "7 Code", 
              "10 Code",  "32 Code")

# Read in table
fileName <- paste(path, "AllGenosmod.csv", sep="")
sampleTable <- read.table(fileName,
                          header = TRUE,
                          sep = ",",
                          stringsAsFactors=FALSE, # Strings like "Cork 4" are not designated as factors
                          check.names=FALSE) # Names left as they are, no dots inserted

# Get the isolates
ourIsolates <- getIsolates(path, "corkisolates.txt")

# Initialise vectors to store info
locationIsolates <- rep(NA, length(ourIsolates))
vntrIsolates <- rep(NA, length(ourIsolates))

# Check if isolate in table, if is, pull out location data and VNTR type
for(row in 1:nrow(sampleTable)){
  for(index in 1:length(ourIsolates)){
    if(ourIsolates[index] %in% sampleTable[row, "MAP Isolates (TB14-00)"] == TRUE){
      locationIsolates[index] <- sampleTable[row, "Herd Location"]
      vntrIsolates[index] <- paste(sampleTable[row, colHeads], collapse = ":")
    }  
  }
}

# Convert VNTR types to INMV identifiers
# Define vector to store new types
inmvIsolates <- rep(NA, length(vntrIsolates))

for(vntr in 1:length(vntrIsolates)){
  if(vntrIsolates[vntr] == "NA:NA:NA:NA:NA:NA:NA:NA" 
     || is.na(vntrIsolates[vntr]) == TRUE){
    inmvIsolates[vntr] <- "NA"
    
  }else if(vntrIsolates[vntr] == "4:2:3:3:2:2:2:8"){
    inmvIsolates[vntr] <- "INMV 1"
  
  }else if(vntrIsolates[vntr] == "3:2:3:3:2:2:2:8"){
    inmvIsolates[vntr] <- "INMV 2"
    
  }else if(vntrIsolates[vntr] == "4:1:3:3:2:2:2:8"){
    inmvIsolates[vntr] <- "INMV 116"
  
  }else if(vntrIsolates[vntr] == "3:2:3:3:2:2:1:8"){
    inmvIsolates[vntr] <- "INMV 3"
  
  }else if(vntrIsolates[vntr] == "4:2:3:3:2:2:3:8" 
           || vntrIsolates[vntr] == "4:2:3:3:2:7:2:8"){
    inmvIsolates[vntr] <- "? INMV 1 ?"
  }
  
}

# Get contaminated isolates
contaminatedIsolates <- getIsolates(path, "contaminateds.txt")

# Get growing isolates
growingIsolates <- getIsolates(path, "growing.txt")

# Get little growth isolates
lgIsolates <- getIsolates(path, "littlegrowing.txt")

# Make a status vector and fill with contaminations if on the list
status <- rep(NA, length(ourIsolates))

for(isolate in 1:length(ourIsolates)){
  if(ourIsolates[isolate] %in% contaminatedIsolates) {
    status[isolate] <- "Contaminated"
  } else if(ourIsolates[isolate] %in% growingIsolates) {
    status[isolate] <- "Growth"
  } else if(ourIsolates[isolate] %in% lgIsolates) {
    status[isolate] <- "Little Growth"
  } else {
    status[isolate] <- "No Growth"
  }
}

# Get rid of duplicate info
ourIsolates[96] <- "5342 3C2"
status[97] <- "NA"

# Create a dataframe to store info
isolateData <- data.frame(ourIsolates, locationIsolates, inmvIsolates, status)

# Write data to table
cat(capture.output(isolateData), file = "isolate_summary_26_03_18.txt", sep = "\n")
write.table(isolateData, file = "isolate_summary_26_03_18.txt", col.names = NA)

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