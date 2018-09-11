
repeatValue <- function(valueToRepeat, nTimesRepeated){
  
  vector <- c()
  
  for(index in 1:nTimesRepeated){
    vector[index] <- valueToRepeat
  }
  
  return(vector)
}
