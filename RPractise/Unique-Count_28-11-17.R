# Define a vector with some elements 
aVector <- c("chicken", "fox", "pigeon", "cat", "chicken", "chicken", "cow",
             "chicken", "cat", "pigeon", "cat", "chicken")

# Function to find unique elements in vector and count how many times
# each shows up

countElements <- function(aVector){
  
  # Initialise list to store element name and its count
  elementList <- list()
  
  # Go thru each element in the vector and add name to list if not present
  # if already there, just add 1 to the count
  for(element in aVector){
    
    if(is.null(elementList[[element]]) == TRUE){
      
      elementList[[element]] <- 1
      
    }
    else{
      
      elementList[[element]] <- elementList[[element]] + 1
    }
    
  }
  return(elementList)
}
  