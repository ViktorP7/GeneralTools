# Load libraries/packages
library(openxlsx)
library(scales)

# Set path variable
path <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/SarahProject/CsqRerun/ComCsq_2020-11-27.tsv"

# Read in table of isolates
csqTable <- read.table(path,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings are not designated as factors
                       check.names=FALSE,
                       quote = "") # Names left as they are, no dots inserted

# Remove the empty bogey column
csqTable <- csqTable[,-length(colnames(csqTable))]

# Get rid of shared consequences between ancestor
relevantCsqs <- ridShared2171Csqs(csqTable)

# 16-2171 column can now be removed as it's empty
relevantCsqs <- relevantCsqs[,-length(colnames(relevantCsqs))]

# Remove any rows full of NAs
onlyRelevantCsqs <- relevantCsqs[-which(is.na(relevantCsqs$`CITP-MAP_csq.vcf`) &
                                      is.na(relevantCsqs$ERR037990_csq.vcf) &
                                      is.na(relevantCsqs$ERR037985_csq.vcf) &
                                      is.na(relevantCsqs$`CIT-MAP_csq.vcf`)), ]

# Get the consequence frequencies including fails
csqFreqsWithFails <- findCsqFreqs(onlyRelevantCsqs, FALSE)

# Get the consequence freqs excluding fails
csqFreqsNoFails <- findCsqFreqs(onlyRelevantCsqs, TRUE)

# Get unique consequences
uniqueCsqs <- findCommonCsqs(onlyRelevantCsqs, 3)

# Get consequences shared across only 2 isolates
twoSharedCsqs <- findCommonCsqs(onlyRelevantCsqs, 2)


# Remove rows that are full of NAs
justUniqueCsqs <- uniqueCsqs[-which(is.na(uniqueCsqs$`CITP-MAP_csq.vcf`) &
                              is.na(uniqueCsqs$ERR037990_csq.vcf) &
                              is.na(uniqueCsqs$ERR037985_csq.vcf) &
                              is.na(uniqueCsqs$`CIT-MAP_csq.vcf`)), ]

justTwoSharedCsqs <- twoSharedCsqs[-which(is.na(twoSharedCsqs$`CITP-MAP_csq.vcf`) &
                                    is.na(twoSharedCsqs$ERR037990_csq.vcf) &
                                    is.na(twoSharedCsqs$ERR037985_csq.vcf) &
                                    is.na(twoSharedCsqs$`CIT-MAP_csq.vcf`)), ]

# Get those that are only unique to CIT and CITP
justCITUniques <- justUniqueCsqs[,-c(3,5)]

# Remove any NA only rows
justCITUniques <- justCITUniques[-which(is.na(justCITUniques$`CIT-MAP_csq.vcf`) &
                                  is.na(justCITUniques$`CITP-MAP_csq.vcf`)),]

# Get those that are shared only by CIT and CITP but not the Bryant stuff
justCITShared <- justTwoSharedCsqs[,-c(3,5)]

# Remove NAs 
justCITShared <- justCITShared[-which(is.na(justCITShared$`CIT-MAP_csq.vcf`) &
                                is.na(justCITShared$`CITP-MAP_csq.vcf`)),]

# Merge the CIT shared and unique frames
mergedCITFrame <- rbind(justCITShared, justCITUniques)

# Find the proportions of consequences for consequences unique to CITs inc. fails
CITFreqsWithFails <- findCsqFreqs(mergedCITFrame, FALSE)
CITFreqsNoFails <- findCsqFreqs(mergedCITFrame, TRUE)

# Convert to proportions
csqPropsWithFails <- convertToProp(csqFreqsWithFails)
csqPropsNoFails <- convertToProp(csqFreqsNoFails)
CITPropsWithFails <- convertToProp(CITFreqsWithFails)
CITPropsNoFails <- convertToProp(CITFreqsNoFails)

# Plot barplot to show overlap of overall consequences and selected consequences
x=barplot(as.numeric(csqPropsNoFails[1:10,4]), col = "red", 
        ylab = "Proportion (%)", names.arg = csqPropsNoFails[1:10,1], las = 2, 
        cex.names = 0.6, main = "Proportions of SNP Consequences for CIT", xaxt = "n")
barplot(as.numeric(CITPropsNoFails[1:10,3]), col = alpha("blue", 0.5),
        add = T, yaxt = "n")

# Add legend
legend("topright", legend = c("Differences from Ancestor 16-2171 (out of 154)",
                           "Differences from Closest Node (out of 19)"),
       text.col = c("Red", alpha("Blue", 0.5)),
       bty = "n", cex = 0.8)

# Add slanted labels
text(x[,1], -0.7, srt = 45, adj= 1, 
     xpd = TRUE, labels = csqPropsNoFails[1:10,1] , cex=0.8)

# Write things to files
write.xlsx(mergedCITFrame, "MergedCITLenient.xlsx")
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

# Function to record frequency of consequences present for each isolate
findCsqFreqs <- function(csqTable, failcheck){
  
  # Create matrix to store results
  resultMatrix <- matrix(data = 0, nrow = 11, ncol = ncol(csqTable), 
                         dimnames = list(c(), c("Type", 
                                           names(csqTable[2:ncol(csqTable)]))))
  
  # Fill in the type column
  resultMatrix[,1] <- c("None", "missense", "frameshift",
                        "inframe_deletion", "synonymous",
                        "inframe_insertion", "stop_gained",
                        "stop_lost&frameshift", "stop_lost",
                        "stop_retained", "Sumtotal")
  
  # Loop thru the input dataframe by column
  for(col in 2:ncol(csqTable)){
    
    # Loop thru the rows of the current column
    for(row in 1:nrow(csqTable)){
      
      # Check if is na or one of the types
      if(is.na(csqTable[row,col]) == TRUE){
        
        next
      
      } else {
        
        # Run a for loop thru the types of consequence
        for(type in 1:nrow(resultMatrix)){
          
          # Check if matches 
          if(grepl(resultMatrix[type,1], csqTable[row, col]) == TRUE){
            
            # Check if failchecker enabled
            if(failcheck == TRUE){
              
              if(grepl("FAIL", csqTable[row, col]) == TRUE){
                
                next
              } else{
                
                resultMatrix[type, col] <- as.numeric(resultMatrix[type, col]) + 1
                
              }
            } else{
              
              resultMatrix[type, col] <- as.numeric(resultMatrix[type, col]) + 1
              
            }
          }
        }
      }
    }
  }
  
  # Tally the results at the end
  for(isolate in 2:ncol(resultMatrix)){
    
    resultMatrix[11,isolate] <- sum(as.numeric(resultMatrix[,isolate]))
  }
  
  return(resultMatrix)
}

# Function to convert to proportions
convertToProp <- function(freqtable){
  
  # Create a result matrix to match the input
  propTable <- freqtable
  
  # Loop thru the columns
  for(col in 2:ncol(freqtable)){
    
    # Loop thru each row
    for(row in 1:nrow(freqtable)){
      
      # Do the conversion
      propTable[row,col] <- round(as.numeric(freqtable[row,col]) / as.numeric(freqtable[11,col]) * 100, digits = 2)
    }
  }
  return(propTable)
}

# Function to get rid of shared consequences between isolates and 16-2171 ancestor
ridShared2171Csqs <- function(csqTable){
  
  # Create a data frame to store results using same col and row names from before
  namesOfRows <- csqTable[,1]
  namesOfCols <- colnames(csqTable)
  filteredFrame <- data.frame(matrix(nrow = length(namesOfRows), 
                                     ncol = length(namesOfCols)))
  names(filteredFrame) <- namesOfCols
  filteredFrame[,1] <- namesOfRows
  
  # Loop thru each row of the ancestral isolate and store if nothing in ancestral
  for(row in 1:length(csqTable$`16-2171_csq.vcf`)){
    
    if(is.na(csqTable$`16-2171_csq.vcf`[row]) == TRUE){
      
      # Add row into the filtered df
      filteredFrame[row,] <- csqTable[row,]
    }
  }
  
  return(filteredFrame)
}
