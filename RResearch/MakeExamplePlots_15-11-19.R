# Extract clades with examples of birth counties 
# Requires other script to work

exampleOne <- extract.clade(subbedTree, node = 239)
# Make EU plot
outputFile <- paste("EG1.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleOne, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = c("deepskyblue3","red","deepskyblue3","deepskyblue3","deepskyblue3","deepskyblue3"),
           align.tip.label = FALSE, type="phylogram", cex = 1)

# Add the SNP scale
add.scale.bar(x=25, y = 2, cex = 1)

dev.off()

exampleTwo <- extract.clade(subbedTree, node = 227)
# Make EU plot
outputFile <- paste("EG2.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleTwo, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = c("deepskyblue3","deepskyblue3"),
           align.tip.label = FALSE, type="phylogram", cex = 1)

# Add the SNP scale
add.scale.bar(x=0.5, y = 1.5, cex = 1)

dev.off()

exampleThree <- extract.clade(subbedTree, node = 198)
# Make EU plot
outputFile <- paste("EG3.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleThree, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = "red",
           align.tip.label = FALSE, type="phylogram", cex = 1)

# Add the SNP scale
add.scale.bar(x=0.5, y = 1.0, cex = 1)

dev.off()

exampleFour <- extract.clade(subbedTree, node = 191)
# Make EU plot
outputFile <- paste("EG4.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleFour, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = c("deepskyblue3","red","deepskyblue3","red","deepskyblue3"),
           align.tip.label = FALSE, type="phylogram", cex = 1)

# Add the SNP scale
add.scale.bar(x=5, y = 1.5, cex = 1)

dev.off()

exampleFive <- extract.clade(subbedTree, node = 168)
# Make EU plot
outputFile <- paste("EG5.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleFive, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = "black",
           align.tip.label = FALSE, type="phylogram", cex = 0.7)

# Add the SNP scale
add.scale.bar(x=10, y = 1.5, cex = 1)

dev.off()

exampleAnimal <- extract.clade(subbedTree, node = 148)
# Make EU plot
outputFile <- paste("Animal.pdf", sep="")
pdf(outputFile, height=5, width=5)

# Plot VNTR tree
plot.phylo(exampleAnimal, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = c("deepskyblue3","red","red","red","deepskyblue3","red","red"),
           align.tip.label = FALSE, type="phylogram", cex = 1)

# Add the SNP scale
add.scale.bar(x=2, y = 1.5, cex = 1)

dev.off()

corkAndEU <- extract.clade(euTree, node = 308)

corkAndEU$tip.label[39] <- "17-5652_Cork_10_1_Cork"

corkColours <- function(tiplabel){
  
  # Create vector to store index values of tips to be dropped
  dropVector <- rep(NA,length(tiplabel))
  
  # Loop thru tip labels and drop as required
  for(index in 1:length(tiplabel)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Cork", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "darkgreen"
    } else if(grepl("Scotland", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "royalblue1"
    } else if(grepl("England", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "red"
    }
  }  
  return(dropVector)  
}

corklabels <- function(tiplabel){
  
  # Create vector to store index values of tips to be dropped
  dropVector <- rep(NA,length(tiplabel))
  
  # Loop thru tip labels and drop as required
  for(index in 1:length(tiplabel)){
    
    # Check if part of a word is in the name and assign a colour
    if(grepl("Cork", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "Cork 10"
    } else if(grepl("Scotland", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "Scotland"
    } else if(grepl("England", tiplabel[index]) == TRUE){
      
      dropVector[index] <- "England"
    }
  }  
  return(dropVector)  
}

corkcols <- corkColours(corkAndEU$tip.label)
corklabs <- corklabels(corkAndEU$tip.label)

corklabs[40] <- "Cork 4"

corkAndEU$tip.label <- corklabs
  
outputFile <- paste("Corker.pdf", sep="")
pdf(outputFile, height=15, width=15)
# Make cork plot
plot.phylo(corkAndEU, edge.width = 2, font = 1, label.offset = 0.2, 
           tip.color = corkcols,
           align.tip.label = FALSE, type="phylogram", cex = 1.5)
add.scale.bar(x=20, y = 1.5, cex = 2)
dev.off()
