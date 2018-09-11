# Initialise S, I and R values, and get total population
s = 50
i = 1
r = 0
n = sum(s,i,r)

# Initialise a time amount, rate of infection, rate of resistance
t = 250
roi = 0.05
ror = 0.005

# run model
rundetSIR(s, i, r, n, t, roi, ror)

# Function to run deterministic SIR simulation
rundetSIR <- function(s, i, r, n, t, roi, ror){
  
  # Create data frame to store values in simulation
  sirData <- data.frame(
    times = c(0:t),
    susceptibles = NA,
    infecteds = NA,
    resistants = NA
  )
  
  # Fill in start values for S, I, R
  sirData[1,"susceptibles"] <- s
  sirData[1,"infecteds"] <- i
  sirData[1,"resistants"] <- r
  
  
  # Increment from 1 until the t value
  for(time in 1:t){
    
    # Calculate change in susceptible population
    changeS <- ((roi*i*s)/n)
    
    # Calculate change in infected population
    changeI <- changeS - ror*i
    
    # Calculate change in resistant population
    changeR <- ror*i
    
    # Apply changes to S, I, R values
    s <- s - changeS
    i <- i + changeI
    r <- r + changeR
    
    # Fill table with S, I and R values corresponding to the time slot
    sirData[time+1,"susceptibles"] <- s
    sirData[time+1,"infecteds"] <- i
    sirData[time+1,"resistants"] <- r
  }
  
  # Plot the data frame to show simulation
  plot(x=sirData$times, y=sirData$susceptibles,
       ylim=c(0, max(sirData[, -1])), type="l", las=1, col="red",
       ylab="N", xlab="Time", main = "SIR Deterministic")
  lines(x=sirData$times, y=sirData$infecteds, col= "green")
  lines(x=sirData$times, y=sirData$resistants, col= "blue")
  legend("left", legend=c("S", "I", "R"),
         text.col=c("red", "green", "blue"), bty="n")
  
}