# Define vector 
someNumbers <- c(1,2,3,4,5,6,7,8,9,12,9,10)

# Function to return max of vector
findMax <- function(someNumbers) {
  currentmaxValue <- someNumbers[1] # set current max = first element

  for(index in 2:length(someNumbers)) { # for every element in the vector
    if(someNumbers[index] > currentmaxValue) { # if they are greater than the 1st num
      currentmaxValue = someNumbers[index] # current max is changed to that
    }
    
  }
  return(currentmaxValue)
}
# Function to return index of vector
findMaxIndex <- function(someNumbers) {
  currentmaxValue <- someNumbers[1]
  maxIndex <- 1 # max index set to 1 to correspond to current max
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] > currentmaxValue) {
      currentmaxValue = someNumbers[index]
      maxIndex <- index # index corresponding to max is set in maxIndex
    }
    
  }
  return(maxIndex)
}

# Write append function
emptyVector <- c()
someElement <- 4

appendValues <- function(emptyVector, someElement) {
  
  emptyVector[length(emptyVector) + 1] <- someElement # elements are added in to the
  # next available space in the empty vector
  
  return(emptyVector)
}
  

# Function to return indices of max values in vector (if more than one)
findMaxIndeces <- function(someNumbers) {
  currentmaxValue <- someNumbers[1]
  
  emptyVector <- c() # empty vector to store the indices of max values
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] > currentmaxValue) {
      currentmaxValue <- someNumbers[index]
      
    }
  }
  
  for(index in 2:length(someNumbers)) {
    if(someNumbers[index] == currentmaxValue) { # fellow max values are placed in
      emptyVector[length(emptyVector) + 1] <- index
    }
  }
  
  return(emptyVector)
}

