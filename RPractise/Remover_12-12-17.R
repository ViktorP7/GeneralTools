# Define a string and a character to remove from it
string <- "chicken"
thing <- "c"

# Define function to remove the character from the string or vector
removeChar <- function(input, thing){
  
  # Empty variable to store output
  output <- NULL
  
  # Split the string up if it is a string and not vector
  if(is.character(input) == TRUE){
    input <- strsplit(input, "")[[1]]
  }
  
  # Loop thru each part of the vector and append it to output if it isn't
  # the object you are looking to remove
  for(index in 1:length(input)){
    if(input[index] != thing){
      output <- append(output, input[index])
    } 
  
  }
  # Paste the string if it's a string
  if(is.character(output) == TRUE){
    output <- paste(output, collapse = "")
  }
  
  return(output)

}