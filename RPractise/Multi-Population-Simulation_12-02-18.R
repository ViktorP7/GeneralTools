# Define an input matrix to store any amount of populations
# Migration rates between populations must be placed in a matrix
# Rows and columns in migration matrix must correspond to amount of populations
# Create parameters matrix to store transitions rates between various sim states
# All matrices used must match each other in terms of disease states and no. pops
parameterMatrix <- matrix(c(0, 0.02, 0, 0, 0, 0, 0.05, 0, 
                            0, 0, 0, 0.03, 0, 0, 0, 0), nrow = 4, byrow = TRUE,
                          dimnames = list(c("S", "E", "I", "R"), 
                                          c("S", "E", "I", "R")))
populationMatrix <- matrix(c(2000, 100, 50, 0, 2000, 10, 5, 0, 
                             2000, 20, 10, 0), nrow = 3, byrow = TRUE,
                           dimnames = list(c("pop1", "pop2", "pop3"), 
                                           c("S", "E", "I", "R")))
migrationMatrix <- matrix(nrow = 3, ncol = 3, 
                          c(0, 0.02, 0.01, 0.025, 0, 0.03, 0.06, 0.1, 0), 
                          dimnames = list(c("pop1", "pop2", "pop3"), 
                                          c("pop1", "pop2", "pop3")))
timer <- c(0, 500)


# Function to run stochastic SEIR simulation for 2 populations
runNPopStocSim <- function(parameterMatrix, populationMatrix, migrationMatrix,
                              timer){
  
  # Initialise time variable
  time <- timer[1]
  
  # Create a vector to store the sum of the pop matrix
  sumVector <- c()
  
  # Get the sums of each population and append them to vector
  for(row in 1:nrow(populationMatrix)){
    
    sumVector <- append(sumVector, sum(populationMatrix[row,]))
  }
  
  # Extract names from the popmatrix to put into state matrix
  namePop <- row.names(populationMatrix)
  nameCondition <- colnames(populationMatrix)
 
  # New vector to store combined names to put into the state matrix 
  comboNames <- c()
  
  # Loop thru names and combine to form new ones
  for(popname in namePop){
    
    for(condname in nameCondition){
      
      currentname <- paste(popname, ",", condname, sep = "")
      
      comboNames <- append(comboNames, currentname)
    }
  }
  
  # Initial matrix needs to be rearranged and placed into the state matrix
  # Initialise the state matrix, must be pops*conditions for amount of cols
  stateMatrix <- matrix(ncol = (ncol(populationMatrix) * nrow(populationMatrix) +1)
                        , dimnames = list(c(), c("time", comboNames)))
  stateMatrix[,1] <- time
  
  # Loop thru rows of pop matrix and fill in the values into state
  # Set up a counter for placing into the state
  counter <- 0
  
  for(row in 1:nrow(populationMatrix)){
    
    for(value in populationMatrix[row,]){
      
      counter <- counter + 1
      
      stateMatrix[1, counter + 1] <- value
      
    }
  }
  
  # Need to find locations of "I" columns in the state matrix
  locatI <- grep("I", comboNames) + 1
  
  # Start a counter to keep track of the while loop
  counterWhile <- 0
  
  # Status for breaking or not
  goStatus <- TRUE

  # While loop to go thru the process
  # Conditions - timer must not go over limit, S column in state must be all over 0
  while(time < timer[2] && stateMatrix["I"]>0){
    
    # Update counter
    counterWhile <- counterWhile + 1
    
    # Check if any I's are below 0
    for(index in locatI){
      
      if(stateMatrix[counterWhile, index] > 0){
      
      } else{
        
        # change go status
        
        goStatus <- FALSE
      }
      
    }
    
    # Double check the go status
    if(goStatus == FALSE){
      
      break
    }
    
    # Calculate probabilities for current state
    currProbs <- getProbs(parameters, state1, state2, oneToTwo, twoToOne)
    
    # Get tau using current probability value to find out when next step happens
    tau <- rexp(n=1, rate=sum(currProbs))
    
    # Update the time
    time <- time + tau
    
    # Find out which process happens after tau by doing a sample roll
    samRoll <- sample(length(currProbs),1,prob=currProbs)
    
  }
  
  
  
}

# Define function for probabilities of infection or resistance occuring
getProbs <- function(parameterMatrix, stateMatrix, migrationMatrix, 
                     counterWhile, sumVector, namePop, comboNames){
  
  # List to store locations
  indexList <- list()
  
  # using the condition names vector, grep out indices of these in the state
  for(pop in namePop){
    
    indexPop <- grep(pop, comboNames) + 1
    
    indexList[[pop]] <- indexPop
                             
  }
  
  # Create new list to store the values
  valueList[namePop] <- list()
  
  # Take the indexList and use the indices to pull out values from the state
  for(pop in 1:length(indexList)){
    
    
    # Inner loop to use indices to pull out values from state, place in list 
    for(index in indexList[[pop]]){
      
      # Variable to store values
      currvals <- NULL
      
      currvals <- stateMatrix[counterWhile, index]
      
      valueList[[pop]] <- append(valueList[[pop]],currvals) 
      
      
    }
  }
 
  # Initialise a vector to store the probabilities
  pVector <- c()
  
 
  return(pVector)
}

# Function to calculate the exposure
pExposure <- function(parameterMatrix, stateValueI, stateValueS, pVector){
  
  pVector <- append(pVector, (parameterMatrix["S", "E"] * stateValueI * stateValueS/sum()))
  
  
}