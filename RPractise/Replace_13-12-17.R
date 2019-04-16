# Define a string and a character to remove from it, and its replacement
string <- "chicken"
ridOut <- "c"
placeIn <-  "k"


# Define function to remove the character from the string or vector
replaceChar <- function(input, ridOut, placeIn){
  
  # Split the string up 
  input <- strsplit(input, "")[[1]]
  
  # Initialise empty vector of same length as input
  output <- rep(NA, length(input))
  
  # Loop thru each part of the string and append it to output if it isn't
  # the object you are looking to remove
  for(index in 1:length(input)){
    if(input[index] != ridOut){
      output[index] <- input[index]
    } 
    
  }
  
  # Loop thru the output and put the desired thing where NAs exist
  for(index in 1:length(output)){
    if(is.na(output[index]) == TRUE){
      output[index] <- placeIn
    }
  }
  
  # Paste the string back together
  output <- paste(output, collapse = "")
  
  
  return(output)
  
}