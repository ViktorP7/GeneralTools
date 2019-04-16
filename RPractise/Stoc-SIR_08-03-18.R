# Initialise S, I and R values, and the parameters
parameters <- c(roi=0.05, ror=0.005)
initial <- c(S=50, I=1, R=0)
timer <- c(0, 250)


# run model
runstocSIR(parameters, initialSIR, time)

# Function to run stochastic SIR simulation
runstocSIR <- function(parameters, initialSIR, time){
  
  # Initialise state and time variables
  state <- initial
  time <- timer[1]
  
  # Create data frame to store values in simulation
  sirData <- data.frame(
    times=time,
    susceptibles = state["S"],
    infecteds = state["I"],
    resistants = state["R"]
  )
  
  # Define changes for S, I and R for each process
  processes <- matrix(0, nrow=2, ncol=3, 
                      dimnames = list(c("infection", "resistance"),
                                      c("dS", "dI", "dR")))
  processes["infection","dS"] <- -1
  processes["infection","dI"] <- +1
  processes["resistance", "dI"] <- -1
  processes["resistance", "dR"] <- +1
  
  # Define function for probabilities of infection or resistance occuring
  getProbs <- function(parameters, state){
    pInfection <- (parameters[["roi"]]*state[["S"]]*state[["I"]])/sum(state)
    pResistance <- (parameters[["ror"]]*state[["I"]])
    pVector <- c(pInfection, pResistance)
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
    
    # The roll will either be a 1(infection) or a 2(resistance)
    if(samRoll == 1){
      
      state <- state + processes["infection",]
      
    } else if(samRoll == 2){
      
      state <- state + processes["resistance",]
    }
    
    # Store values into the output data frame
    sirData <- rbind(sirData, c(time, state))
  }
  
  # Get max time for plot limit
  maxTime <- max(sirData["times"])
  
  
  plot(x=sirData$times, y=sirData$susceptibles,
       ylim=c(0, max(sirData[, -1])), type="l", las=1, col="red",
       ylab="N", xlab="Time", main = "SIR Stochastic")
  lines(x=sirData$times, y=sirData$infecteds, col= "green")
  lines(x=sirData$times, y=sirData$resistants, col= "blue")
  legend("left", legend=c("S", "I", "R"),
         text.col=c("red", "green", "blue"), bty="n")
  
  
}