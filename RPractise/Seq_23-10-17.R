# Define a start value
startValue <- 0

# Define an end value
endValue <- 30

# Define a by value
byValue <- 5

genSeq <- function(startValue, endValue, byValue){

  # Calculate the length of the vector to be used
  fillerInitial <- (endValue - startValue)/byValue
  fillerLength <- floor(fillerInitial) +  1

  # Initialise empty vector for function to fill
  fillerVector <- rep(NA, fillerLength)

  # Place the start value into position one of vector
  fillerVector[1] <- startValue

  # Add values in increments of by to vector
  for(index in 2:length(fillerVector)) {
    Increment <- fillerVector[index - 1] 
    fillerVector[index] <- Increment + byValue
  }
  return(fillerVector)
}