# Define a vector
someVector <- c(2,3,4,5,6,4,3,4,5,7,8,7,5,5,6)

medianFunction <- function(someVector){
  # Sort the vector 
  sortedVector <- sort(someVector, decreasing=FALSE)

  # Find the median position
  medianPosition <- length(sortedVector)/2
  if(medianPosition %% 2 == 0) {
    realMed <- (medianPosition + (medianPosition + 1)) / 2
  } else {
    realMed <- ceiling(medianPosition)
  }
  med <- sortedVector[realMed]
  return(med)
}