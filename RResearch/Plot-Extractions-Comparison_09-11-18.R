# Get necessary packages
library(dplyr)
library(scales)

# Set path and read in the file
path <- "C:/Users/UCD/Documents/Lab/Nanodrop/extraction-comparison_9-11-18.csv"

extractionTable <- read.table(path, 
                          header=TRUE, 
                          sep=",", 
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)

# Set the concentration threshold
concentrationThreshold <- 200

# Replace minus values with zeroes by looping thru rows
extractionTable[extractionTable$`Total DNA ng` < 0, "Total DNA ng"] <- 0

# Plot out DNA vs OD
plotDNAvODinBeadsvPhenolvCTAB(extractionTable, concentrationThreshold)

# Plot 260/280 ratio vs 260/230
plotQualityinBeadsvPhenolVCTAB(extractionTable)

# Filter out outlier ratios
filteredTable <- extractionTable %>%
  filter(`260/280` <= 2.5 & `260/280` >= 1.5 & `260/230` <= 2.5 & `260/230` >= 1.5 & `Total DNA ng` >= 200)

# Plot 260/280 ratio vs 260/230 excluding outliers
plotQualityinBeadsvPhenolVCTAB(filteredTable)

# Change the x and y axis scale
axis(side=1, at=seq(1.5,2.5,by=0.1))
axis(side=2, at=seq(1.5,2.5,by=0.1), las=1)

# Plot boxplots of 260/280 that passed criteria 200ng
plotQualityBoxplots(filteredTable, "260/280")

# Plot boxplots of 260/230 that passed criteria 200ng
plotQualityBoxplots(filteredTable, "260/230")

###############
###FUNCTIONS###
###############
# Function to plot DNA vs OD
plotDNAvODinBeadsvPhenolvCTAB <- function(extractionTable, concentrationThreshold){
  
  
  # Plot OD600 reading against total DNA concentration
  plot(x=extractionTable[, "OD600nm"], y=extractionTable[, "Total DNA ng"], # Provide X and Y coordinates
       main="DNA Extracted from 500ul liquid", # Set title
       xlab="OD600nm", # Label X axis
       ylab="Total DNA (ng)", # Label y axis
       log = "y", # make y axis logarithmic
       xaxt = "n",
       las=1, # Change angle of tick labels to be horizontal
       bty="n", # Remove box around plot
       col=ifelse(grepl("ctab", extractionTable[,"Sample ID"]) , rgb(0,1,0, 1),
                  ifelse(grepl("beads", extractionTable[,"Sample ID"]), rgb(1,0,0,1), rgb(0,0,1,1))),
       pch=ifelse(grepl("ctab", extractionTable[,"Sample ID"]) , 0,
                  ifelse(grepl("beads", extractionTable[,"Sample ID"]), 1, 2))) 
  
  # Change the x axis scale
  axis(side=1, at=seq(0,2,by=0.1), las = 2)
  
  # Add a legend
  legend("bottomright", # Set the position of the legend
         legend=c("Beads", "Phenol", "CTAB"), # Set the names of the legend
         text.col=c("red", "blue", "green"), # Set the colours of the legend names
         pch=c(1,2,0), # Set the point shape
         col=c("red", "blue","green"), # Set the point colours
         bty="n") # Remove the box around
  
  # Add in DNA concentration threshold
  abline(h=concentrationThreshold, 
         col="grey", # Set the colour of the line
         lty=2) # Make the line dashed
}

# Function to plot quality ratios
plotQualityinBeadsvPhenolVCTAB <- function(extractionTable){
  plot(x=extractionTable[, "260/280"], y=extractionTable[, "260/230"], # Provide X and Y coordinates
       main="DNA Quality", # Set title
       xlab="Absorbance ratio at 260/280nm", # Label X axis
       ylab="Absorbance ratio at 260/230nm", # Label y axis
       las=1, # Change angle of tick labels to be horizontal
       bty="n", # Remove box around plot
       col=ifelse(grepl("ctab", extractionTable[,"Sample ID"]) , rgb(0,1,0, 1),
                  ifelse(grepl("beads", extractionTable[,"Sample ID"]), rgb(1,0,0,1), rgb(0,0,1,1))),
       pch=ifelse(grepl("ctab", extractionTable[,"Sample ID"]) , 0,
                  ifelse(grepl("beads", extractionTable[,"Sample ID"]), 1, 2)))
  
  # Add a legend
  legend("bottomleft", # Set the position of the legend
         legend=c("Beads", "Phenol", "CTAB"), # Set the names of the legend
         text.col=c("red", "blue", "green"), # Set the colours of the legend names
         pch=c(1,2,0), # Set the point shape
         col=c("red", "blue","green"), # Set the point colours
         bty="n") # Remove the box around  
  
  
  # Add in lines for inner thresholds
  abline(v=c(1.8, 2.0), 
         col="black", # Set the colour of the line
         lty=2) # Make the line dashed
  
  # Add in lines for outer thresholds
  #abline(h=c(1.85,2.25), v=c(1.75, 2.05), 
  #       col="grey", # Set the colour of the line
  #       lty=2) # Make the line dashed
  
  # Draw a box to highlight desired quality
  #rect(xleft = 1.8, ybottom = 2, xright = 2, ytop = 2.2)
  
  # Draw a box to highlight acceptable bounds for 260/230
  #rect(xleft = 1.75, ybottom = 1.85, xright = 2.05, ytop = 2.25, border = "grey")
  
}

# Function to plot quality as boxplots
plotQualityBoxplots <- function(extractionTable, ratio){
  
  # Do a boxplot
  boxplot(extractionTable[grep("ctab", extractionTable[,"Sample ID"]),ratio], 
          extractionTable[grep("beads", extractionTable[,"Sample ID"]),ratio],
          extractionTable[grep("phenol", extractionTable[,"Sample ID"]),ratio],
          names = c("CTAB", "Beads", "Phenol"),
          main = paste(ratio, "quality ratios for DNA methods"),
          ylab = paste(ratio, "absorbance ratio"))
  
  stripchart(extractionTable[grep("ctab", extractionTable[,"Sample ID"]),ratio],
             add = TRUE, at = 1, col = alpha("green", 0.5), 
             pch = 0, vertical = TRUE, method = "jitter")
  stripchart(extractionTable[grep("beads", extractionTable[,"Sample ID"]),ratio],
             add = TRUE, at = 2, col = alpha("red", 0.5), 
             pch = 1, vertical = TRUE, method = "jitter")
  stripchart(extractionTable[grep("phenol", extractionTable[,"Sample ID"]),ratio],
             add = TRUE, at = 3, col = alpha("blue", 0.5), 
             pch = 2, vertical = TRUE, method = "jitter")
}