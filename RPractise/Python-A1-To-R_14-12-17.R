# Designate path 
path <- "path/to/files"

# Designate the file name
fileName <- "inputA.txt"

# Call function to read in fasta file
fastaFile <- readFasta(path, fileName)

# Call function to get the names of the fasta sequences
seqNames <- getFastaNames(fastaFile)

# Call function to get the sequences
seqs <- getFastaSeqs(fastaFile)

# Call function to get lengths of all sequences
seqLenghts <- getSeqLengths(seqs)



###FUNCTIONS###

# Function to read in the input line by line
readFasta <- function(path, fileName){
  
  # Paste together the path and file name
  input <- paste(path, fileName, sep = "")
  
  # Read in the fasta file line by line
  fastaInput <- readLines(input, warn = FALSE)
  
  return(fastaInput)
}

# Function to store names of sequences
getFastaNames <- function(fastaFile){
  
  # Initialise vector to store the names
  fastaNames <- c()
  
  # Loop thru each line and check if it has a ">" denoting name
  for(line in fastaFile){
    
    # Variable to store current line
    currline <- NULL
    
    # Split each line up
    currline <- strsplit(line, "")[[1]]
    
    # Check for ">"
    if(currline[1] == ">"){
      
      fastaNames <- append(fastaNames, line)
    }
  }
  return(fastaNames)
}

# Function to store sequences themselves
getFastaSeqs <- function(fastaFile){
  
  # Initialise vector to store the seqs
  fastaSeqs <- c()
  
  # Loop thru each line and check if it has a ">" denoting name
  for(line in fastaFile){
    
    # Variable to store current line
    currline <- NULL
    
    # Split each line up
    currline <- strsplit(line, "")[[1]]
    
    # Check for ">"
    if(currline[1] == ">"){
      
      next
      
    }else{
      
      fastaSeqs <- append(fastaSeqs, line)
    }
  }
  return(fastaSeqs)
}

# Function to calculate the length of each sequence, ignoring spaces
getSeqLengths <- function(seqs){
  
  # Initialise vector to store the seq lengths
  lengthsVector <- NULL
  
  # Loop thru the sequences
  for(sequence in seqs){
    
    # Variable to store each split sequence
    currseq <- NULL
    
    # Split up each sequence
    currseq <- strsplit(sequence, "")[[1]]
    
    # Grep out the spaces from the sequence
    currseq <- grep("-", currseq, invert = TRUE)
    
    # Calculate the length of currseq and append it to the lengths vector
    lengthsVector <- append(lengthsVector, length(currseq))
  }
  
  return(lengthsVector)
}