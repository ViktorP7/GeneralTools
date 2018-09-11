# Load libraries/packages
library(xlsx)

# Set path variable
path <- "path/to/tsv"

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
uniqueCsqs <- findCommonCsqs(csqTable, 3)

# Get consequences shared across only 2 isolates
twoSharedCsqs <- findCommonCsqs(csqTable, 2)

# Get consequences shared across 3 isolates
threeSharedCsqs <- findCommonCsqs(csqTable, 1)

# Get consequences common to all
sharedCsqs <- findCommonCsqs(csqTable, 0)

# Remove rows that are full of NAs
justUniqueCsqs <- uniqueCsqs[-which(is.na(uniqueCsqs$`CITP-MAP_csq.vcf`) &
                              is.na(uniqueCsqs$ERR037990_csq.vcf) &
                              is.na(uniqueCsqs$ERR037985_csq.vcf) &
                              is.na(uniqueCsqs$`CIT-MAP_csq.vcf`)), ]

justTwoSharedCsqs <- twoSharedCsqs[-which(is.na(twoSharedCsqs$`CITP-MAP_csq.vcf`) &
                                    is.na(twoSharedCsqs$ERR037990_csq.vcf) &
                                    is.na(twoSharedCsqs$ERR037985_csq.vcf) &
                                    is.na(twoSharedCsqs$`CIT-MAP_csq.vcf`)), ]

justThreeSharedCsqs <- threeSharedCsqs[-which(is.na(threeSharedCsqs$`CITP-MAP_csq.vcf`) &
                                      is.na(threeSharedCsqs$ERR037990_csq.vcf) &
                                      is.na(threeSharedCsqs$ERR037985_csq.vcf) &
                                      is.na(threeSharedCsqs$`CIT-MAP_csq.vcf`)), ]

allSharedCsqs <- sharedCsqs[-which(is.na(sharedCsqs$`CITP-MAP_csq.vcf`) &
                              is.na(sharedCsqs$ERR037990_csq.vcf) &
                              is.na(sharedCsqs$ERR037985_csq.vcf) &
                              is.na(sharedCsqs$`CIT-MAP_csq.vcf`)), ]

# Get those that are only unique to CIT and CITP
justCITUniques <- justUniqueCsqs[,-c(3,4)]

# Remove any NA only rows
justCITUniques <- justCITUniques[-which(is.na(justCITUniques$`CIT-MAP_csq.vcf`) &
                                  is.na(justCITUniques$`CITP-MAP_csq.vcf`)),]

# Get those that are shared only by CIT and CITP but not the Bryant stuff
justCITShared <- justTwoSharedCsqs[,-c(3,4)]

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
