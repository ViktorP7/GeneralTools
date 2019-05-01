# Load the timevis package for making a GANTT chart
library(timevis)
library(htmlwidgets)

# Read in the GANTT table
gantt <- read.table("C:/Users/UCD/Documents/R Tools/PhD_GANTT_26-03-19.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Set the date formats
gantt$start <- as.Date(gantt$start, format="%d/%m/%Y")
gantt$end <- as.Date(gantt$end, format="%d/%m/%Y")

# Remove NA rows
naRows <- which(is.na(gantt$id))
if(length(naRows) > 0){
  gantt <- gantt[-naRows, ]
}

# Note the groupings in the gannt chart
#groups <- data.frame(
#  id = c("General", "Chapter 4", "Chapter 5", "Chapter 6"),
#  content = c("General", "Chapter 4", "Chapter 5", "Chapter 6")
#)

# Set the options for timevis
options <- list(editable = FALSE, align = "center", multiselect=TRUE)

# Plot the GANTT chart - Click Export -> "Save as Web Page..." to save as interactive html!
timevis(gantt, options=options)

# Save the timeline to file
saveWidget(timevis(gantt, options=options), 
           "gantt1-05-19.html", selfcontained=FALSE)

