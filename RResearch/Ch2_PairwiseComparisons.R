library(phytools)
library(ape)
library(ggplot2)
library(egg)

correctnames <- c("13-2711","13-4284","13-4516","13-4573","13-4959","13-5083","13-6823","13-8781","14-3230",
                  "14-5172","14-5298","14-6194","14-6342","14-6538","14-6677","14-7626","14-7739","14-8576",
                  "14-9252", "16-5154", "Ref")

assoc <- cbind(c("13-2711","13-4284","13-4516","13-4573","13-4959","13-5083","13-6823","13-8781","14-3230",
                 "14-5172","14-5298","14-6194","14-6342","14-6538","14-6677","14-7626","14-7739","14-8576",
                 "14-9252", "16-5154", "Ref"),
               c("13-2711","13-4284","13-4516","13-4573","13-4959","13-5083","13-6823","13-8781","14-3230",
                 "14-5172","14-5298","14-6194","14-6342","14-6538","14-6677","14-7626","14-7739","14-8576",
                 "14-9252", "16-5154", "Ref"))



# Unfiltered trees
sim.nf <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.orig_15-09")
jcp.nf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.jNF_14-09')
snpgenie.nf <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.sNF_14-09")
snippy.nf <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.snippyNF_14-06")

# Change some tiplabels
sim.nf$tip.label[18] <- "Ref"
snpgenie.nf$tip.label[5] <- "Ref"
sim.nf$tip.label <- fixthenames(correctnames, sim.nf$tip.label)
jcp.nf$tip.label <- fixthenames(correctnames, jcp.nf$tip.label)
snpgenie.nf$tip.label <- fixthenames(correctnames, snpgenie.nf$tip.label)
snippy.nf$tip.label <- fixthenames(correctnames, snippy.nf$tip.label)


unfilt1 <- cophylo(sim.nf,jcp.nf,assoc,rotate=TRUE)
plot(unfilt1, link.lwd=2,link.col="red",fsize=0.6, scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "UNFILTERED; Sim vs JCP")

unfilt2 <- cophylo(sim.nf,snpgenie.nf,assoc,rotate=TRUE)
plot(unfilt2, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "UNFILTERED; Sim vs SNiPgenie")

unfilt3 <- cophylo(sim.nf,snippy.nf,assoc,rotate=TRUE)
plot(unfilt3, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "UNFILTERED; Sim vs snippy")

# Filtered and filtered default trees
sim.df <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.origSF04-11")
jcp.df <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.JCPDefault_14-09')
snpgenie.df <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.variants")
snippy.df <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.snippyDefault_14-06")

jcp.sf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.JCSFNP20-09')
snpgenie.sf <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.SGSF20-09")
snippy.sf <- read.tree(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.snippySF_14-09")

jcp.psf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.JCSFP20-09')
jcp.rsf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/RAxML_bipartitions.JCSFRNP20-09')

# Change some tiplabels
sim.df$tip.label[17] <- "Ref"
snpgenie.df$tip.label[5] <- "Ref"
sim.df$tip.label <- fixthenames(correctnames, sim.df$tip.label)
jcp.df$tip.label <- fixthenames(correctnames, jcp.df$tip.label)
snpgenie.df$tip.label <- fixthenames(correctnames, snpgenie.df$tip.label)
snippy.df$tip.label <- fixthenames(correctnames, snippy.df$tip.label)

snippy.sf$tip.label <- fixthenames(correctnames, snippy.sf$tip.label)
snpgenie.sf$tip.label[15] <- "Ref"
snpgenie.sf$tip.label <- fixthenames(correctnames, snpgenie.sf$tip.label)
jcp.sf$tip.label <- fixthenames(correctnames, jcp.sf$tip.label)
jcp.psf$tip.label <- fixthenames(correctnames, jcp.psf$tip.label)
jcp.rsf$tip.label <- fixthenames(correctnames, jcp.rsf$tip.label)

filtx <- cophylo(sim.df,jcp.df,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filtx, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Default; Sim vs JCP")

filty <- cophylo(sim.df,snpgenie.df,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filty, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Default; Sim vs SNiPgenie")

filtz <- cophylo(sim.df,snippy.df,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filtz, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Default; Sim vs snippy")

filt1 <- cophylo(sim.df,snippy.sf,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt1, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Standard; Sim vs snippy")

filt2 <- cophylo(sim.df,snpgenie.sf,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt2, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Standard; Sim vs SNiPgenie")

filt3 <- cophylo(sim.df,jcp.sf,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt3, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Standard; Sim vs JCP")

filt6 <- cophylo(snippy.nf,snippy.df, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt6, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Snippy; Unfiltered vs Default")

filt7 <- cophylo(snpgenie.nf,snpgenie.df, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt7, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "SNiPGenie; Unfiltered vs Default")

filt8 <- cophylo(jcp.nf, jcp.df, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt8, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "JCP; Unfiltered vs Default")


filt9 <- cophylo(jcp.sf, jcp.psf, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt8, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "JCP; No Prox vs Prox")

filt10 <- cophylo(jcp.sf, jcp.rsf, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt8, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "JCP; No Rescue vs Rescue")

filt15 <- cophylo(sim.nf,sim.df, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt15, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.2,0.2))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
text(x=0,y=-0.1, "Sim: Unfiltered vs Filtered")


# Function to correct pipeline names
fixthenames <- function(correct, query){
   
  outvector <- query
  
  for(position in correct){
    
    place <- grep(position, query)
    
    outvector[place] <- position
  }
  return(outvector)
}