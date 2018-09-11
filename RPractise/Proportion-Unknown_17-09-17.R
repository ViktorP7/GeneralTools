# Define sequence list
# As this function focuses on getting the proportion of "N"s in a sequence, the lengths
# of the sequences don't necessarily have to be equal
sequences <- list(
  seqA = "ATGCNNTGC",
  seqB = "ACNGNTGCATNT",
  seqC = "NNATCGC"
)

# Function to calculate proportion of unknown nucleotides in sequences
calculatePropotionUnknown <- function(sequences) {

  # Get an array of keys for the list
  # Keys provide access to the elements of the list, in this case, the names of
  # the sequences act as the keys to the relevant sequences
  keys <- names(sequences)
  
  # Convert sequences to character vectors
  # This for loop is able to access each sequence, first thru to whatever the last
  # might be, as defined by the length (or amount) of keys available. Thus, a task
  # can be performed upon in sequence. Index is used to call up each key, which in
  # turn calls up each sequence. The seqs, which prior to this particular code
  # appear as strings, are now broken into individual characters by the strsplit
  # "" is used as part of the split in order to split each character with a space.
  for(index in 1:length(keys)) {
    
    sequences[[keys[index]]] <- strsplit(sequences[[keys[index]]], split="")[[1]]
  }
  
  # Initialise a vector to store the results of the "N" proportion for each sequence
  proportionsN <- rep(NA, length(keys))
  
  # Count the number of "N" values in each sequence
  # This can be done using a for loop that cycles through each sequence, and adds a 1
  # to the countN variable for each "N" detected 
  # However, a nested for loop is then needed to calculate proportion "N" for each seq
  # and store this value in the proportionsN vector. 
  for(sequenceIndex in 1:length(keys)) {
#   This is the short way of achieving the same function    
#   proportionsN[sequenceIndex] <- length(which(sequences[[keys[sequenceIndex]]] == "N")) /
#      length(sequences[[keys[sequenceIndex]]])
    
    # Get the current sequence
    # This is a shortcut for accessing each sequence in the list
    currentSequence <- sequences[[keys[sequenceIndex]]]
    
    # Initialise a variable to store the "N" count for each sequence
    # Must be initialised inside the for loop as it is reset for each sequence
    countN <- 0
    
    # nucleotideIndex accesses each individual nucleotide in each sequence, and thus
    # each is examined for if it is N or not. 
    for(nucleotideIndex in 1:length(currentSequence)){
      
      if(currentSequence[nucleotideIndex] == "N") {
        
        countN <- countN + 1
      }
    }
    
    # Divide by length of sequence to get proportion of "N"
    # Note that nucleotideIndex does not exist outside the nested for loop, and thus 
    # cannot be used anymore in an upper level of code
    proportionsN[sequenceIndex] <- countN / length(currentSequence)
    
  }
  
  return(proportionsN)
}