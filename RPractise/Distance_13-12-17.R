# Define 4 variables corresponding to x1, y1, x2, y2
xOne <- 6
yOne <- 7
xTwo <- 9
yTwo <- 2

# Function to return euclidean distance
distPoints <- function(xOne, yOne, xTwo, yTwo){
  
  # Create variable to store distance
  distVal <- NULL
  
  # Calculate and store the distance
  distVal <- sqrt((xOne - xTwo)^2 + (yOne - yTwo)^2)
  
  return(distVal)
  
}