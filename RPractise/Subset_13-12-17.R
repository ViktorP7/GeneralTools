# Define vector, start and end positions
numbers <- c(17,22,33,45,56,61)
start <- 2
end <- 5


# Function to return subset of the vector
getsubset <- function(numbers, start, end){
  
  # Define vector to store subset
  subvector <- c()
  
  # Loop thru input using start and end as indices
  for(index in start:end){
    subvector <- append(subvector, numbers[index])
  }
  
  return(subvector)
}