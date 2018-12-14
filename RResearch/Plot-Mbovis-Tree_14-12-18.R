# Load packages
library(ape)
library(geiger)
library(phytools)

# Set path variables
pathTree <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/Bovis-FASTQs/vcfFiles/BovisTree.tree"


# Read in tree
TheTree <- read.tree(pathTree)

# Plot the tree
plot.phylo(TheTree, edge.width = 0.2, font = 1, label.offset = 0.01)
nodelabels()
tiplabels()

pdf("BovisTree.pdf", width=20, height=20)

# Root the tree at reference 
rootree <- root(TheTree, "Ref-1997")

# Convert branch lengths to SNP values
rootree$edge.length <- rootree$edge.length * 390

# Get the floored values o the lengths
flooredSNPs <- floor(rootree$edge.length)

# Drop the ref
droppedTree <- drop.tip(rootree, "Ref-1997")

# Plot the new tree
plot.phylo(droppedTree, edge.width = 0.2, font = 1, label.offset = 0.001, 
           align.tip.label = TRUE, type="phylogram", cex = 0.9)

# Add the SNP scale
add.scale.bar(cex = 1.0)

# Add SNP lengths to each edge
edgelabels(text = flooredSNPs, adj = c(0.5, -0.25), frame = "none", cex = 0.5)

dev.off()


