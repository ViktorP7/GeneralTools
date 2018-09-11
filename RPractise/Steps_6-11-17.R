n <- 10

for(i in 1:n){
  Sys.sleep(0.5)
  
  cat(paste("\rStep ", i, " of ", n , sep=""))
  
}