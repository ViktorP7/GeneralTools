# Define sequence string
someSeq <- "ATGCCGTCA"

# Function to get complement of a sequence
getComplement <- function(someSeq) {

  # Convert sequence to character vector
  splitSeq <- strsplit(someSeq, split = "")[[1]]
  
  # Define complementary bases
  bases <- c("A", "T", "G", "C", "N")
  names(bases) <- c("T", "A", "C", "G", "N") # dictionary allows conversion
  
  # bases <- list(
  #   "A" = "T",
  # ...
  #)
  # bases[[]]
  
  # Initialise variable to store converted sequence
  complementSeq <- rep(NA,length(splitSeq))
  
  # For loop to swap out the bases in input sequence
  for(index in 1:length(splitSeq)) { # for each base 
    complementSeq <- bases[splitSeq] # convert each base to complement, store in vector
  }
  
  
  return(paste(complementSeq, collapse ="")) # turns individual bases to string again

}