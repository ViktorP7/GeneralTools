# Define vector 
inputVector <- c(2,3,4,4,3,4)

meanFunction <- function(inputVector){
  
  # Initialise a variable to store the calculated mean
  mean <- NULL

  # Check that the input is numeric
  if(is.numeric(inputVector) == FALSE){
    print("Input is not in the correct format")
  }else{
  
    # Initialise a variable to store the sum of the input vector
    sum <- 0
  
    # Calculate the mean
    for(index in 1:length(inputVector)){
      sum <- inputVector[index] + sum
    }

    mean <- sum / length(inputVector)
  }
  
  return(mean)
}
