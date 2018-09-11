# Read in data from file
table <- read.table("path/to/csv",
                       header = TRUE,
                       sep = ",")
# Tell where file is, write path
# Comma separated, has column headers

# Function to plot difference between before and after values
beforeVafter <- function(table) {
  
  if(ncol > 3 || ncol < 3){
    print("Table has incorrect amount of columns")
  } else {
  

  # Plot data from table

  plot(x=table[, "Before"], y=table[, "After"], # input data
       xlim = c(0,1), ylim = c(0,1), # Setting the axis limits
       main = "Difference between before and after", # Changing title
      xlab = "Before", ylab = "After", # Setting axis names
      las=1) # Rotating y axis values to horizontal position

  lines(x=c(0,1), y=c(0,1), col = "red", lty = 2)
  # add line going diagonal y=x, line of equality, function lines
  # make line red and dashed
  
  # For each row in the table, plot line from point where x=before, y=after to 
  # line of equality where x=y=before, make line green and dashed
  for(row in 1:nrow(table)) {
    
    lines(x=c(table[row, "Before"], table[row, "Before"]), 
          y=c(table[row, "After"], table[row, "Before"]),
          col = "green", lty = 2 )
  
    Sys.sleep(0.05)
    }
  }
}

#while loop way of doing it
#row <- 0
#while(row <= nrow(table)) {
#  
#  row <- row + 1
#  
#  lines(x=c(table[row, "Before"], table[row, "Before"]),
#        y=c(table[row, "After"], table[row, "Before"]),
#        col = "green", lty = 2 )
#  
#}


