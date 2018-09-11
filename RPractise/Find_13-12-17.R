# Define vector and value to find in vector
somenums <- c(1,2,3,4,5,43,2,2,1,2,44,3,21,1)
target <- 2

# Function to find locations of target in a vector
findIndeces <- function(somenums, target) {
  
  # Empty vector to store indices
  emptyVector <- c() 
  
  # Loop thru vector and see if any values match the target, store index in
  # vector if they do
  for(index in 1:length(somenums)) {
    if(somenums[index] == target) { 
      emptyVector[length(emptyVector) + 1] <- index
    }
  }
  
  return(emptyVector)
}