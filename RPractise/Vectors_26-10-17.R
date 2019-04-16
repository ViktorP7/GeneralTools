# Define two vectors
vectorA <- c(1,2,3,4,5,6,7) 
vectorB <- c(2,3,4,6,7,4,3)

vectorAddition <- function(vectorA, vectorB) {
  
  # Check that vector lengths are equal for both vectors
  if(length(vectorA) != length(vectorB)) {
    print("Vectors not equal length, cannot carry out function")
  } else {
    
    # Initialise an output vector
    vectorC <- rep(NA, length(vectorA))
    
    # Add elements of one vector to the other
    for(index in 1:length(vectorA)) {
      
      
      vectorC[index] <- vectorA[index] + vectorB[index]
    }
    return(vectorC)
  }
}

vectorSubtraction <- function(vectorA, vectorB) {
  
  # Check that vector lengths are equal for both vectors
  if(length(vectorA) != length(vectorB)) {
    print("Vectors not equal length, cannot carry out function")
  } else {
    # Add elements of one vector to the other
    for(index in 1:length(vectorA)) {
      vectorC[index] <- vectorA[index] - vectorB[index]
    }
    return(vectorC)
  }
}

vectorMulti <- function(vectorA, vectorB) {
  
  # Check that vector lengths are equal for both vectors
  if(length(vectorA) != length(vectorB)) {
    print("Vectors not equal length, cannot carry out function")
  } else {
    # Add elements of one vector to the other
    for(index in 1:length(vectorA)) {
      vectorC[index] <- vectorA[index] * vectorB[index]
    }
    return(vectorC)
  }
}

vectorDiv <- function(vectorA, vectorB) {
  
  # Check that vector lengths are equal for both vectors
  if(length(vectorA) != length(vectorB)) {
    print("Vectors not equal length, cannot carry out function")
  } else {
    # Add elements of one vector to the other
    for(index in 1:length(vectorA)) {
      vectorC[index] <- vectorA[index] / vectorB[index]
    }
    return(vectorC)
  }
}
