# Define a vector and a value
vectorA <- c(1,2,3)
value <- 1

inVector <- function(vector, value){
  # Variable to record presence
  present <- FALSE

  # Cycle through the vector and check if the value is present
  for(index in 1:length(vector)){
    if(vector[index] == value){
      present <- TRUE
      break
    }
  }

  return(present)
}  
