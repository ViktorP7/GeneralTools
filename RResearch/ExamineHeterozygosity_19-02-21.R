# Load libraries/packages
library(openxlsx)
library(scales)

# Set path variable
pathAD <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/MAPAncestors/vcfFiles/AD_2021-03-15.tsv"

# Read in table of isolates
adTable <- read.table(pathAD,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors=FALSE, # Strings are not designated as factors
                       check.names=FALSE,
                       quote = "") # Names left as they are, no dots inserted

# Remove the empty bogey column
adTable <- adTable[,-length(colnames(adTable))]

# Save plots in pdf
pdf("HetSummary_15-03-21.pdf")

layout(matrix(c(1:length(colnames(adTable))-1)), nrow=4, ncol=4, byrow=TRUE)

# Loop thru the isolates and create a plot for each
for(index in 2:length(colnames(adTable))){
  
  valvector <- adTable[,index]
  
  valvector <- valvector[!is.na(valvector)]
  
  totalpos <- length(valvector)
  
  plot(x=NA, y=NA,
       xlim = c(0, totalpos), 
       ylim = c(0, 1),
       main = colnames(adTable)[index], xlab = "Position", ylab = "Proportion of Supporting Reads",
       bty = "n")
  
  counter = 0

  
  for(position in 1:length(valvector)){
    
    curRef <- as.numeric(strsplit(valvector[position], split = ",")[[1]][1])
    curAlt <- as.numeric(strsplit(valvector[position], split = ",")[[1]][2])
    totalDP <- sum(curRef, curAlt)
    propRef <- round(curRef/totalDP, digits = 2)
    propAlt <- round(curAlt/totalDP, digits = 2)
    
    if(propAlt > 0 && propRef > 0){
      
      counter = counter + 1
    }
    
    points(x=position, y=propRef, pch = 20, col = alpha("blue", 0.5))
    points(x=position, y=propAlt, pch = 20, col = alpha("red", 0.5))
                         
  }
  
  perHet <- round((counter/totalpos)*100, digits = 2)
  legend("right", legend = c("Alt","Ref"), bty="n", cex=1.5,
         pch = c(20,20), col = c(alpha("red", 0.9), alpha("blue",0.9)),
         title = paste(perHet, "% places are het"))
  
}

dev.off()

# Investigate the most commonly shared het positions
commonShared <- data.frame()

for(row in 1:nrow(adTable)){
  
  if(sum(is.na(adTable[row,])) > 100){
    
    next
  }else if(sum(grepl("0,", adTable[row,])) > 50){
    
    next
  }else{
    
    commonShared <- rbind(commonShared, adTable[row,])
  }
}