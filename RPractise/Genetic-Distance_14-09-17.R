# Define sequence strings
seqA <- "ACTCCGTG"
seqB <- "ACTCGCTG"

# Function to calculate genetic difference between two sequences
calculateGenDiff <- function(seqA, seqB) {
  
  # Convert sequences to character vectors
  seq1 <- strsplit(seqA, split="")[[1]]
  seq2 <- strsplit(seqB, split="")[[1]]

    # Initialise variable to record number of differences
    genDiff <- 0
  
    # Check sequences are the same length
    if(length(seq1) != length(seq2)) {
      print("Sequences are not of equal length")
    } else {
      
      # Count number of genetic differences
      for(index in 1:length(seq1)) {
        if(seq1[index] != seq2[index]) {
          genDiff <- genDiff + 1
        }
      }
    }
  return(genDiff)  
}    
