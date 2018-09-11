# Initialise S, E, I and R values, and the parameters
parameters <- c(roi=0.01, ror=0.0005, inc=0.005)
initial <- c(S=1500, E=20, I=10, R=0)
timer <- c(0, 2000)



# Function to run stochastic SEIR simulation
runstocSEIR <- function(parameters, initialSEIR, time){

  # Initialise state and time variables
  state <- initial
  time <- timer[1]
  
  # Create data frame to store values in simulation
  seirData <- data.frame(
    times=time,
    susceptibles = state["S"],
    exposeds = state["E"],
    infecteds = state["I"],
    resistants = state["R"]
  )
  
  # Define changes for S, E, I and R for each process
  processes <- matrix(0, nrow=3, ncol=4, 
                      dimnames = list(c("exposure","infection", "resistance"),
                                      c("dS", "dE", "dI", "dR")))
  processes["exposure", "dS"] <- -1
  processes["exposure","dE"] <- +1
  processes["infection","dE"] <- -1
  processes["infection","dI"] <- +1
  processes["resistance", "dI"] <- -1
  processes["resistance", "dR"] <- +1
  
  # Define function for probabilities of infection or resistance occuring
  getProbs <- function(parameters, state){
    pExposure <- (parameters[["roi"]]*state[["S"]]*state[["I"]])/sum(state)
    pInfection <- (parameters[["inc"]]*state[["E"]])
    pResistance <- (parameters[["ror"]]*state[["I"]])
    pVector <- c(pExposure,pInfection, pResistance)
    return(pVector)
  }
  
  
  # While loop to go thru the process
  while(time < timer[2] && state["I"]>0){
    
    # Calculate probabilities for current state
    currProbs <- getProbs(parameters, state)
    
    # Get tau using current probability value to find out when next step happens
    tau <- rexp(n=1, rate=sum(currProbs))
    
    # Update the time
    time <- time + tau
    
    # Find out which process happens after tau by doing a sample roll
    samRoll <- sample(length(currProbs),1,prob=currProbs)
    
    # The roll will either be a 1(exposure), 2(infection) or a 3(resistance)
    if(samRoll == 1){
      
      state <- state + processes["exposure",]
      
    } else if(samRoll == 2){
      
      state <- state + processes["infection",]
    
    } else if(samRoll == 3){
      
      state <- state + processes["resistance",]
    }
    
    # Store values into the output data frame
    seirData <- rbind(seirData, c(time, state))
  }
  
  # Get max time for plot limit
  maxTime <- max(seirData["times"])
  
  
  plot(x=seirData$times, y=seirData$susceptibles,
       ylim=c(0, max(seirData[, -1])), type="l", las=1, col="red",
       ylab="N", xlab="Time", main = "SEIR Stochastic")
  lines(x=seirData$times, y=seirData$exposeds, col= "purple")
  lines(x=seirData$times, y=seirData$infecteds, col= "green")
  lines(x=seirData$times, y=seirData$resistants, col= "blue")
  legend("left", legend=c("S", "E", "I", "R"),
         text.col=c("red", "purple", "green", "blue"), bty="n")
  
  
}