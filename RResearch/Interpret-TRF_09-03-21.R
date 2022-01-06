library(Biostrings)
library(DECIPHER)

location = "C:/Users/UCD/Desktop/UbuntuSharedFolder/ResultsTRF/"
pathMeta <- "C:/Users/UCD/Documents/Lab/CVRL MAP/ForComparisonMetaSep2021Format.csv"

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

finalFrame <- data.frame(matrix(nrow = length(files), ncol = 1))

colnames(finalFrame) <- c("INMV Code")

rownames(finalFrame) <- sapply(strsplit(files, "\\."), "[", 1)

for(row in 1:nrow(resultFrame)){
  
  finalFrame[row,1] <- paste(resultFrame[row,],collapse = "")
}


write.csv(finalFrame, file = paste(format(Sys.time(), "%Y-%m-%d"),"VNTRresults.csv",sep = "_"))

metaTable <- read.table(pathMeta,
                           header = TRUE,
                           sep = ",",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)

# Check real vs computed VNTR
comparison <- compareVNTRs(resultFrame, metaTable, 196)
  

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
  
  if(reNum == 8){
  
    alifilter <- 100
    indfilter <- 20
  } else if(reNum == 5) {
    
    alifilter <- 90
    indfilter <- 3
  } else if(reNum == 6){
    
    alifilter <- 90
    indfilter <- 10
  } else{
    alifilter <- 90
    indfilter <- 5
  }
  for(row in 1:nrow(inTable)){
    
  # Filter by repeat length, alignment score and indel presence between repeats
    if(inTable[row,5] == reLength[reNum] && inTable[row,8] >= alifilter && inTable[row,7] <= indfilter){
    
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
      } else if(length(bestmatch) >=2 && reNum ==8){
        
        topscorealn <- which(max(interTable[bestmatch,8]) == interTable[bestmatch, 8])
        
        outTable[currFile, reNum] <- rounder(interTable[bestmatch[topscorealn],4])
      } else if(length(bestmatch) >=2){
        
        topscoreadj <- which(max(interTable[bestmatch,6]) == interTable[bestmatch, 6])
        topscoreind <- which(min(interTable[bestmatch,7]) == interTable[bestmatch, 7])
        
        if(length(topscoreadj) == 1 && length(topscoreind) == 1 && topscoreadj == topscoreind){
          
          toptop <- topscoreadj
        } else if(topscoreadj %in% topscoreind){
          toptop = topscoreadj
        }else{
          toptop = topscoreind[1]
        }
        
        outTable[currFile, reNum] <- rounder(interTable[bestmatch[toptop],4])
      }
    } else if(length(bestmatch) ==1){
      
      outTable[currFile, reNum] <- rounder(interTable[bestmatch[1],4] + interTable[nearestmatch[1],4])
        
    } else{
        
      topscoreadj <- which(max(interTable[bestmatch,6]) == interTable[bestmatch, 6])
      topscoreind <- which(min(interTable[bestmatch,7]) == interTable[bestmatch, 7])
      
      if(length(topscoreadj) == 1 && length(topscoreind) == 1 && topscoreadj == topscoreind){
        
        toptop <- topscoreadj
      } else if(topscoreadj %in% topscoreind){
        toptop = topscoreadj
      }else{
        toptop = topscoreind[1]
      }
      
      outTable[currFile, reNum] <- rounder(interTable[bestmatch[toptop],4])
    }
  }
  
  return(outTable)
}

# Convert results to types and compare
compareVNTRs <- function(resTable, meta, qLen){
  
  outVect <- rep(NA, qLen)

  for(row in 1:qLen){
    
    if(meta[row,14] == 1){
      
      currvec <- c(4,2,3,3,2,2,2,8)
    } else if(meta[row,14] == 2){
      
      currvec <- c(3,2,3,3,2,2,2,8)
      
    }else if(meta[row,14] == 3){
      
      currvec <- c(3,2,3,3,2,2,1,8)
      
    } else if(meta[row,14] == 13){
      
      currvec <- c(2,2,3,3,2,2,2,8)
      
    }else if(meta[row,14] == 116){
      
      currvec <- c(4,1,3,3,2,2,2,8)
      
    }
    
    diffVec = resTable[row,] - currvec
    
    outVect[row] = length(which(diffVec != 0))
    
  }
  return(outVect)
}
