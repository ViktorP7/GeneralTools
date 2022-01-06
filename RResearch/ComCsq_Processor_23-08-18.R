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
relevantCsqs <- csqTable[,-length(colnames(csqTable))]


# Check if any rows full of NAs
rowstorem = which(rowSums(is.na(relevantCsqs)) == (ncol(relevantCsqs)-1))

# Get the consequence frequencies including fails
csqFreqsWithFails <- findCsqFreqs(relevantCsqs, FALSE)

# Get the consequence freqs excluding fails
csqFreqsNoFails <- findCsqFreqs(relevantCsqs, TRUE)

# Get unique consequences
uniqueCsqs <- findCommonCsqs(relevantCsqs, (ncol(relevantCsqs)-2))

# Get consequences shared across only 2 isolates
twoSharedCsqs <- findCommonCsqs(relevantCsqs, (ncol(relevantCsqs)-3))


# Remove rows that are full of NAs
urowstorem = which(rowSums(is.na(uniqueCsqs)) == (ncol(uniqueCsqs)-1))
justUniqueCsqs <- uniqueCsqs[-urowstorem, ]

trowstorem = which(rowSums(is.na(twoSharedCsqs)) == (ncol(twoSharedCsqs)-1))
justTwoSharedCsqs <- twoSharedCsqs[-trowstorem, ]

# Get those that are only unique to CIT and CITP
justCITUniques <- justUniqueCsqs[,-c(2,5,6,7)]

# Remove any NA only rows
justCITUniques <- justCITUniques[-which(is.na(justCITUniques$`CIT-MAP_csq.vcf`) &
                                  is.na(justCITUniques$`CITP-MAP_csq.vcf`)),]

# Get those that are shared only by CIT and CITP but not the Bryant stuff
justCITShared <- justTwoSharedCsqs[,-c(2,5,6,7)]

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
        cex.names = 0.6, main = "Proportions of Consequences for CIT", xaxt = "n")


# Add legend
legend("topright", legend = "Differences from K10 (out of 333)",
       text.col = "red",
       bty = "n", cex = 0.8)

# Add slanted labels
text(x[,1], -0.7, srt = 45, adj= 1, 
     xpd = TRUE, labels = csqPropsNoFails[1:10,1] , cex=0.8)

x=barplot(as.numeric(CITPropsWithFails[1:10,3]), col = alpha("blue", 0.8),ylab = "Proportion (%)",names.arg = CITPropsWithFails[1:10,1]
        , las=2, cex.names = 0.6, main = "Proportions of Consequences for CIT", xaxt = "n")

# Add legend
legend("topright", legend = "Differences from Nearest Node (out of 20*)",
       text.col = "blue",
       bty = "n", cex = 0.8)

# Add slanted labels
text(x[,1], -0.7, srt = 45, adj= 1, 
     xpd = TRUE, labels = CITPropsWithFails[1:10,1] , cex=0.8)

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

