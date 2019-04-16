# Load libraries/packages
library(xlsx)

# Set path variable
path <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/ExtendedConsequences/ComCsq_2018-12-17.tsv"

# Read in table of isolates
csqTable <- read.table(path,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings are not designated as factors
                       check.names=FALSE) # Names left as they are, no dots inserted

# Remove the empty bogey column
csqTable <- csqTable[,-length(colnames(csqTable))]

# Create variable to store how many NAs to look for 
sumNAs <- 3

# Get unique consequences
uniqueCsqs <- findCommonCsqs(csqTable, 11)

# Get consequences shared across only 2 isolates
twoSharedCsqs <- findCommonCsqs(csqTable, 10)

# Get consequences shared across 3 isolates
threeSharedCsqs <- findCommonCsqs(csqTable, 9)

# Get consequences common to all
sharedCsqs <- findCommonCsqs(csqTable, 0)

# Remove rows that are full of NAs
justUniqueCsqs <- removeNAs(uniqueCsqs)
justTwoSharedCsqs <- removeNAs(twoSharedCsqs)
justThreeSharedCsqs <- removeNAs(threeSharedCsqs)


# Get those that are only unique to CIT and CITP
justCITUniques <- justUniqueCsqs[,c(1,7,9)]

# Remove any NA only rows
justCITUniques <- justCITUniques[-which(is.na(justCITUniques$`CIT-MAP_csq.vcf`) &
                                          is.na(justCITUniques$`CITP-MAP_csq.vcf`)),]

# Get those that are shared only by CIT and CITP but not the Bryant stuff
justCITShared <- justTwoSharedCsqs[,c(1,7,9)]

# Remove NAs 
justCITShared <- justCITShared[-which(is.na(justCITShared$`CIT-MAP_csq.vcf`) &
                                        is.na(justCITShared$`CITP-MAP_csq.vcf`)),]

# Write things to files
write.xlsx(justCITUniques, "CITUniquesLenient.xlsx")
write.xlsx(justCITShared, "CITSharedLenient.xlsx")
#############
##FUNCTIONS##
#############

# Function to find shared/unshared consequences
findCommonCsqs <- function(csqTable, sumNAs){
  
  # Create a data frame to store results using same col and row names from before
  namesOfRows <- csqTable[,1]
  namesOfCols <- colnames(csqTable)
  filteredFrame <- data.frame(matrix(nrow = length(namesOfRows), 
                                     ncol = length(namesOfCols)))
  names(filteredFrame) <- namesOfCols
  filteredFrame[,1] <- namesOfRows
  
  # Loop thru each row of the table
  for(row in 1:nrow(csqTable)){
    
    # Initiate an NA counter
    naCount <- 0
    
    # Loop thru each item and row and update counter if it is an NA
    for(item in csqTable[row,]){
      
      if(is.na(item)){
        
        naCount <- naCount + 1
        
      }
    }
    
    # Check if the amount of NAs equals the specified sum
    if(naCount == sumNAs){
      
      # Add row into the filtered df
      filteredFrame[row,] <- csqTable[row,]
      
    }
  }
  
  return(filteredFrame)
}

# Function to remove NA rows
removeNAs <- function(csqTable){
  
  # Create a data frame to store results using same col and row names from before
  namesOfCols <- colnames(csqTable)
  filteredFrame <- data.frame(matrix(,ncol = length(namesOfCols)))
  names(filteredFrame) <- namesOfCols
  
  # Loop thru each row of the table
  for(row in 1:nrow(csqTable)){
    
    # Initiate an NA counter
    naCount <- 0
    
    # Loop thru each item and row and update counter if it is an NA
    for(item in csqTable[row,]){
      
      if(is.na(item)){
        
        naCount <- naCount + 1
        
      }
    }
    
    # Check if the amount of NAs equals the specified sum
    if(naCount < (length(namesOfCols)-1)){
      
      # Add row into the filtered df
      filteredFrame <- rbind(filteredFrame, csqTable[row,])
      
    }
  }
  
  filteredFrame <- filteredFrame[-1,]
  
  return(filteredFrame)
}
