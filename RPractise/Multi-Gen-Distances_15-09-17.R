# Define sequence strings
# Sequences are put into a list, allowing the use of keys, and enabling an x 
# amount of sequences to be processed by the function, as opposed to just a few
# The sequences start off as strings "", and are separated by commas in the list
# Throughout the code, elements of a list as accessed by [[]]
# However, elements of a vector or matrix only need a single square brackets

sequences <- list(
  seqA = "ACTCCNTG",
  seqB = "ACTNGCTN",
  seqC = "TATCGGAC"
)

#sequences <- c("ACTCCGTG", "ACTCGCTG", "TATCGGAC")
#names <- c("seqA", "seqB", "seqC")

# Function to calculate genetic difference between sequences
calculateMultiGenDist <- function( sequences) {
  
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
  
  # Initialise a vector to store lengths of sequences
  # NA is used as the values of the vector are yet unknown, but the length of the
  # vector is defined as the length of the keys, in other words the amount of seqs
  sequenceLengths <- rep(NA, length(keys))
  
  # Get the length of each sequence
  # This for loop has access to each key, and thus the index of each key is used
  # to access each respective sequence, calculating the length. The empty vector 
  # coded previously is now filled with the length of each sequence.
  for(index in 1:length(keys)) {
    sequenceLengths[index] <- length(sequences[[keys[index]]])
  }
  
  # Initialise a matrix to store the genetic distances between the input sequences
  #     1   2   3
  # 1  1,1 1,2 1,3
  # 2  2,1 2,2 2,3
  # 3  3,1 3,2 3,3
  # The rows and columns in the matrix are defined by the length of the keys,
  # which indicate the number of sequences. A matrix is needed because several
  # sequences are being compared to one another.
  genDistances <- matrix(nrow = length(keys), ncol = length(keys))
  
  # Check the length of sequences is the same
  # An if statement is used to check that all sequences are of the same length by
  # measuring the amount of unique sequence lengths. If this amount is not equal 
  # to 1, the function is terminated and an error message is shown. Else, the 
  # function carries on.
  if(length(unique(sequenceLengths)) != 1) {
    print("Sequences are not of equal length") 
  }else{
    
    # Get each sequence
    # The structure below is a nested for loop. This means that for each cycle
    # of the upper for loop, the lower full loop goes thru all cycles, ie
    # (1(1,2,3)2(1,2,3)3(1,2,3), etc), allowing each sequence to be compared to
    # all others, including itself. For the first loop, index is defined as "I"
    # so as to avoid conflict with the second loop, where index is "J". In both 
    # cases, index is the length of keys, ie amount of sequences. For each loop,
    # the sequences being compared are assigned seqI and seqJ respectively.
    for(indexI in 1:length(keys)) {
      
      seqI <- sequences[[keys[indexI]]]
      
      # Compare the sequence you got to every other sequence
      for(indexJ in 1:length(keys)) {
        
        seqJ <- sequences[[keys[indexJ]]]
        
        # Skip when indexI = indexJ
        # Since the matrix produced at the end of this will be symmetric, it is
        # easier to leave half the values blank, as they will be the same either
        # way. THis if statement makes the matrix fill only when indexI (rows) is 
        # smaller than indexJ (columns). This saves computing power in large tasks

        if(indexI < indexJ) {
          
          # Calculate genetic distances between sequences
          # The function created previously, whereby two sequences are compared
          # to get the genetic distance between them, is used to calculate the
          # genetic distance between seqs I and J. The values are filled into
          # relevant matrix positions, denoted by [indexI, indexJ] in genDistances
          # [indexI, indexJ] is equaled to [indexJ, indexI] to fill the spaces
          # previously left blank. As this is just a copying of values, less 
          # time is devoted to the calculations, making the function more efficient
          genDistances[indexI,indexJ] <- calculateGenDiff(seqI, seqJ)
          genDistances[indexJ,indexI] <- genDistances[indexI, indexJ]
          # The genDistances matrix can now be returned to show values
        } 
      }  
    }
  }
  return(genDistances)
}  
  
  
  
calculateGenDiff <- function( seqA, seqB) {
  
  # Assumes sequences are the same length
  # Sequences must be arrays of characters
  
  # Initialise variable to record number of differences
  genDiff <- 0
    
  # Count number of genetic differences
  for(index in 1:length(seqA)) {
    if(seqA[index] != "N" && seqB[index] != "N" && seqA[index] != seqB[index]) {
      genDiff <- genDiff + 1
    }
  }
    
  return(genDiff)  
}   
