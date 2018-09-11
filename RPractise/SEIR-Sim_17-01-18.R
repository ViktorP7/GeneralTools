# Initialise S, E, I and R values, and get total population
s = 250
e = 20
i = 10
r = 0
n = sum(s,e,i,r)

# Initialise a time amount, rate of infection, rate of resistance, incubation
t = 2000
roi = 0.01
ror = 0.0005
inc = 0.005

# Function to run deterministic SEIR simulation
rundetSEIR <- function(s, e, i, r, n, t, roi, ror, inc){
  
  # Create data frame to store values in simulation
  sirData <- data.frame(
    times = c(0:t),
    susceptibles = NA,
    exposeds = NA,
    infecteds = NA,
    resistants = NA
  )
  
  # Fill in start values for S, E, I, R
  sirData[1,"susceptibles"] <- s
  sirData[1,"exposeds"] <- e
  sirData[1,"infecteds"] <- i
  sirData[1,"resistants"] <- r

  
  # Increment from 1 until the t value
  for(time in 1:t){
    
    # Calculate change in susceptible population
    changeS <- ((roi*i*s)/n)
    
    # Calculate change in exposed population
    changeE <- changeS - inc*e
    
    # Calculate change in infected population
    changeI <- inc*e - ror*i
    
    # Calculate change in resistant population
    changeR <- ror*i
    
    # Apply changes to S, E, I, R values
    s <- s - changeS
    e <- e + changeE
    i <- i + changeI
    r <- r + changeR
    
    # Fill table with S, E, I and R values corresponding to the time slot
    sirData[time+1,"susceptibles"] <- s
    sirData[time+1, "exposeds"] <- e
    sirData[time+1,"infecteds"] <- i
    sirData[time+1,"resistants"] <- r
  }

  # Plot the data frame to show simulation
  plot(x=sirData$times, y=sirData$susceptibles,
       ylim=c(0, max(sirData[, -1])), type="l", las=1, col="red",
       ylab="N", xlab="Time", main = "SEIR Deterministic")
  lines(x=sirData$times, y=sirData$exposeds, col="purple")
  lines(x=sirData$times, y=sirData$infecteds, col= "green")
  lines(x=sirData$times, y=sirData$resistants, col= "blue")
  legend("left", legend=c("S", "E", "I", "R"),
         text.col=c("red", "purple", "green", "blue"), bty="n")
 
}
