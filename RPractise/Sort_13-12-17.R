# Define a vector
vector <- c(1,2,4353,6,7,65,5,756,45,43,67,67,54,654,77,99,1,3543,6)

# Function to sort the vector
sortIt <- function(vector){
  
  cat("Starting ...\t\t\t\t\t\t\t\t\t\t")
  # Bring the max element to the end of vector using outer loop
  for(index in 1:(length(vector)-1)){
    
    # Nested loop to do comparison
    for(nIndex in 1:(length(vector)-index)){
      
      Sys.sleep(0.5)
      cat(paste("\r", paste(vector, collapse=",")))
      
      # If current element is bigger than next, swap them over
      if(vector[nIndex] > vector[nIndex + 1]){
        
        current <- vector[nIndex]
        vector[nIndex] <- vector[nIndex + 1]
        vector[nIndex + 1] <- current
      }
    }
  }
  cat("\n")
  
  return(vector)
  
}