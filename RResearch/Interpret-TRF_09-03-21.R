library(Biostrings)
library(DECIPHER)

location = "C:/Users/UCD/Desktop/UbuntuSharedFolder/ResultsTRF/"

# Script to parse TRF output and write results to file
files <- list.files(path = location, pattern = "\\.dat$")

# Initialise dataframe to store result
resultFrame <- data.frame(matrix(nrow = length(files), ncol = 8))

colnames(resultFrame) <- c("292", "X3", "25", "47", "3", "7", "10", "32")

rownames(resultFrame) <- sapply(strsplit(files, "\\."), "[", 1)

repLengths <- c(53,53,58,35,27,22,55,18)

repSeqs <- c("GTCATCTGAGCCGCTCCTCCTCATCGCTGCGCTCTGCATCGTCGTCGGCGCGA", "CGCGCGTATGACGATGCAGAGCGTAGCGATGAGGAGGAATCGCGCCGATGACC",
             "AGCCCGGCGACGAAGCATTCCGCGCAGCGGGATGCGCAGGAGCGGGCGATTGGGGTGG",
             "GAGACGGGCGCCCTGATGAGACCCGCGCCCGG", "GCTGGCCGCGATGCCGACGCTGCCGGA", "CGAAATATTCGCCGTGAGAACA",
             "CCGCGGCGATCGCAAGCGCGGCGGAGCCGGGCGCCGCGGGTCGCCGCCATCGACA",
             "CCGCCGCCGGGCTACGGC")

# Read in the files
for(index in 1:length(files)){
  
  # Read in table of isolates
  daTable <- read.table(paste(location,files[index], sep = ""),
                        header = FALSE,
                        sep = " ",
                        stringsAsFactors=FALSE, # Strings are not designated as factors
                        check.names=FALSE,
                        quote = "") # Names left as they are, no dots inserted


  for(num in 1:8){
    resultFrame <- repeatScanner(repLengths, repSeqs, daTable, resultFrame, num, index)
  } 
}

write.csv(resultFrame, file = paste(format(Sys.time(), "%Y-%m-%d"),"VNTRresults.csv",sep = "_"))
  

# Function that rounds up at .5
rnd <- function(x) trunc(x+sign(x)*0.5)

# Choose a rounder
rounder <- function(val){
  
  if(val == 2.5)
    
    outval <- rnd(val)
  else{
    
    outval <- round(val)
  }
  return(outval)
}
  
# Function to scan for repeats
repeatScanner <- function(reLength, reSeq, inTable, outTable, reNum, currFile){
  
  interTable <- data.frame()
  
  for(row in 1:nrow(inTable)){
  
    if(inTable[row,5] == reLength[reNum] && inTable[row,8] >= 100){
    
      interTable <- rbind(interTable, inTable[row,])
    }
  }
  
  if(nrow(interTable) == 1){
    
    outTable[currFile, reNum] <- rounder(interTable[1,4])
  } else {
      
    scores <- rep(NA,nrow(interTable))
      
    for(spot in 1:nrow(interTable)){
        
      aligno <- AlignSeqs(myXStringSet = c(DNAStringSet(reSeq[reNum]), DNAStringSet(interTable[spot,15])), verbose = FALSE)
        
      alignmat <- DistanceMatrix(aligno)
        
      scores[spot] <- alignmat[2]
    }
      
    bestmatch <- which(scores == min(scores))
      
    nearestmatch <- which(scores <= min(scores)+0.01)
      
    nearestmatch <- nearestmatch[-which(nearestmatch == bestmatch)]
      
    if(length(nearestmatch) == 0){
      
      if(length(bestmatch) == 1){
      
        outTable[currFile, reNum] <- rounder(interTable[bestmatch,4])
        
      } else if(length(bestmatch) >= 2 && reNum == 1) {
        
        outTable[currFile, reNum] <- rounder(interTable[bestmatch[1],4] + interTable[bestmatch[2],4])
      }
    } else if(length(bestmatch) ==1){
      
      outTable[currFile, reNum] <- rounder(interTable[bestmatch[1],4] + interTable[nearestmatch[1],4])
        
    } else{
        
      outTable[currFile, reNum] <- rounder(interTable[bestmatch[1],4] + interTable[bestmatch[2],4])
    }
  }
  
  return(outTable)
}
