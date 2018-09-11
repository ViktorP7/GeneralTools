# Define a vector
somevector <- c(1,2,3,4,5)

# Function to reverse the vector
vectorReversal <- function(somevector){
  
  # Create an empty vector to store the reverse 
  reversedVector <- c()
  
  # Go thru each element of the vector in reverse and fill it into the rev
  for(index in length(somevector):1){
    reversedVector <- append(reversedVector, somevector[index])
  }
  
  return(reversedVector)
}