# Initialise S, E, I and R values for populations 1 and 2, and the parameters
# Also initialise the migration rates between the two populations
parameters <- c(roi=0.05, ror=0.005, inc=0.01)
initial1 <- c(S=1200, E=10, I=10, R=0)
initial2 <- c(S=1000, E=0, I=5, R=100)
timer <- c(0, 500)
oneToTwo <- 0.006
twoToOne <- 0.003


# Function to run stochastic SEIR simulation for 2 populations
runTwoPopStocSEIR <- function(parameters, initial1, initial2, time, 
                              oneToTwo, twoToOne){
  
  # Initialise state and time variables
  state1 <- initial1
  state2 <- initial2
  time <- timer[1]
  
  
  # Create data frame to store values in simulation
  seirData <- data.frame(
    times=time,
    susceptibles1 = state1["S"],
    exposeds1 = state1["E"],
    infecteds1 = state1["I"],
    resistants1 = state1["R"],
    susceptibles2 = state2["S"],
    exposeds2 = state2["E"],
    infecteds2 = state2["I"],
    resistants2 = state2["R"]
  )
  
  # Define changes for S, E, I and R for each process
  processes <- matrix(0, nrow=5, ncol=4, 
                      dimnames = list(c("exposure","infection", "resistance",
                                        "migrationSin", "migrationSout"),
                                      c("dS", "dE", "dI", "dR")))
  processes["exposure", "dS"] <- -1
  processes["exposure","dE"] <- +1
  processes["infection","dE"] <- -1
  processes["infection","dI"] <- +1
  processes["resistance", "dI"] <- -1
  processes["resistance", "dR"] <- +1
  processes["migrationSin", "dS"] <- +1
  processes["migrationSout", "dS"] <- -1
  
  # While loop to go thru the process
  while(time < timer[2] && state1["I"]>0 && state2["I"]>0){
    
    # Calculate probabilities for current state
    currProbs <- getProbs(parameters, state1, state2, oneToTwo, twoToOne)
    
    # Get tau using current probability value to find out when next step happens
    tau <- rexp(n=1, rate=sum(currProbs))
    
    # Update the time
    time <- time + tau

    # Find out which process happens after tau by doing a sample roll
    samRoll <- sample(length(currProbs),1,prob=currProbs)
    
    # The roll will either be a 1(exposure1), 2(infection1), 3(resistance1),
    # 4(exposure2), 5(infection2), 6(resistance2), 7(migration12) or 8(migration21)
    if(samRoll == 1){
      
      state1 <- state1 + processes["exposure",]
      
    } else if(samRoll == 2){
      
      state1 <- state1 + processes["infection",]
      
    } else if(samRoll == 3){
      
      state1 <- state1 + processes["resistance",]
      
    } else if(samRoll == 4){
        
      state2 <- state2 + processes["exposure",]
      
    } else if(samRoll == 5){
        
      state2 <- state2 + processes["infection",]
        
    } else if(samRoll == 6){
        
      state2 <- state2 + processes["resistance",]
        
    } else if(samRoll == 7){
        
      state1 <- state1 + processes["migrationSout",]
      state2 <- state2 + processes["migrationSin",]
    
    } else if(samRoll == 8){
      
      state1 <- state1 + processes["migrationSin",]
      state2 <- state2 + processes["migrationSout",]
    } 
    
    # Store values into the output data frame
    seirData <- rbind(seirData, c(time, state1, state2))
  }
  
  
  plot(x=seirData$times, y=seirData$susceptibles1,
       ylim=c(0, max(seirData[, -1])), type="l", las=1, col="red",
       ylab="N", xlab="Time", main = "SEIR Stochastic")
  lines(x=seirData$times, y=seirData$exposeds1, col= "purple")
  lines(x=seirData$times, y=seirData$infecteds1, col= "green")
  lines(x=seirData$times, y=seirData$resistants1, col= "blue")
  lines(x=seirData$times, y=seirData$susceptibles2, col= "pink")
  lines(x=seirData$times, y=seirData$exposeds2, col= "violet")
  lines(x=seirData$times, y=seirData$infecteds2, col= "grey")
  lines(x=seirData$times, y=seirData$resistants2, col= "black")
  legend("left", legend=c("S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2"),
         text.col=c("red", "purple", "green", "blue", "pink", "violet", "grey", "black"), bty="n")
  
  
}

# Define function for probabilities of infection or resistance occuring
getProbs <- function(parameters, stateA, stateB, oneToTwo, twoToOne){
  
  pExposure1 <- (parameters[["roi"]]*stateA[["I"]]*
                   stateA[["S"]])/sum(stateA)
  
  pInfection1 <- parameters[["inc"]]*stateA[["E"]]
  
  pResistance1 <- parameters[["ror"]]*stateA[["I"]]
  
  pExposure2 <- (parameters[["roi"]]*stateB[["I"]]*
                   stateB[["S"]])/sum(stateB)
  
  pInfection2 <- parameters[["inc"]]*stateB[["E"]]
  
  pResistance2 <- parameters[["ror"]]*stateB[["I"]]
  
  pMigrationS12 <- stateA[["S"]]*oneToTwo
  
  pMigrationS21 <- stateB[["S"]]*twoToOne
  
  pVector <- c(pExposure1,pInfection1, pResistance1,
               pExposure1,pInfection1, pResistance1, pMigrationS12, pMigrationS21)
  return(pVector)
}