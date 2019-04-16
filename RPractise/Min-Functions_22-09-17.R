# Define vector 
someNumbers <- c(1,2,3,4,5,6,7,8,9,12,9,10)

# Function to return min of vector
findMin <- function(someNumbers) {
  currentminValue <- someNumbers[1] # set current min = first element
  
  for(index in 2:length(someNumbers)) { # for every element in the vector
    if(someNumbers[index] < currentminValue) { # if they are less than the 1st num
      currentminValue = someNumbers[index] # current min is changed to that
    }
    
  }
  return(currentminValue)
}
# Function to return index of vector
findMinIndex <- function(someNumbers) {
  currentminValue <- someNumbers[1]
  minIndex <- 1 # min index set to 1 to correspond to current min
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] < currentminValue) {
      currentminValue = someNumbers[index]
      minIndex <- index # index corresponding to min is set in minIndex
    }
  }
  return(minIndex)
  
}

# Function to return indices of min values in vector (if more than one)
findMinIndices <- function(someNumbers) {
  currentminValue <- someNumbers[1]
  
  emptyVector <- c() # empty vector to store the indices of min values
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] < currentminValue) {
      currentminValue = someNumbers[index]
    }
  }
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] == currentminValue) { # fellow min values are placed in
      emptyVector[length(emptyVector) + 1] <- index
    }
  }
  
  return(emptyVector)
}

