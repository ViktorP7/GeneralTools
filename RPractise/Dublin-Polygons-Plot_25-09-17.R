# Set path variable
path <- "path/to/polygons" 

# Read in the polygon coordinates
polygonsFile <- paste(path, "PolygonCoords_Dublin.txt", sep="") 
polygonCoords_Dublin <- read.table(polygonsFile, header = TRUE, sep = "\t")

# Calculate the x and y axis limits
rangeX <- range(polygonCoords_Dublin[, "X"])
rangeY <- range(polygonCoords_Dublin[, "Y"])

# Create an empty plot to plot polygon
plot(x=NA, y=NA,
     xlim = rangeX, ylim = rangeY,
     main = "", xlab = "", ylab = "",
     bty = "n", axes = "n")

# Use "polygon" function to add polygon of county onto figure
polygon(x=polygonCoords_Dublin[, "X"],
        y=polygonCoords_Dublin[, "Y"]) # make line a bit thicker

# Place county name roughly in the middle of polygon use "text" function
text(x=mean(polygonCoords_Dublin[, "X"]),
     y=mean(polygonCoords_Dublin[, "Y"]),
     labels = "Dublin")
