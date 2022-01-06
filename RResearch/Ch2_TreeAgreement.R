library(treespace)
library(phytools)
library("adegenet")
library("adegraphics")
library("rgl")
library(ggplot2)
library(rgl)
library(funData)

correctnames <- c("13-2711","13-4284","13-4516","13-4573","13-4959","13-5083","13-6823","13-8781","14-3230",
                  "14-5172","14-5298","14-6194","14-6342","14-6538","14-6677","14-7626","14-7739","14-8576",
                  "14-9252", "16-5154", "Ref")

# Read in all trees
sim_unfilt <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.orig_15-09')
sim_10_filt <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.origSF04-11')
jcp_unfilt <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.jNF_14-09')
jcp_df <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.JCPDefault_14-09')
jcp_sf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.JCSFNP20-09')
sg_uf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.sNF_14-09')
sg_df <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.variants')
sg_sf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.SGSF20-09')
sn_uf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.snippyNF_14-06')
sn_df <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.snippyDefault_14-06')
sn_sf <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bootstrap.snippySF_14-09')

sim_unfilt_b <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.orig_15-09')
sim_10_filt_b <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.origSF04-11')
jcp_unfiltb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.jNF_14-09')
jcp_dfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.JCPDefault_14-09')
jcp_sfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.JCSFNP20-09')
sg_ufb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.sNF_14-09')
sg_dfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.variants')
sg_sfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.SGSF20-09')
sn_ufb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.snippyNF_14-06')
sn_dfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.snippyDefault_14-06')
sn_sfb <- read.tree(file = 'C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/Trees/BS/RAxML_bestTree.snippySF_14-09')

# Fix all their names
sim_unfilt <- runfixnames(100, sim_unfilt, correctnames)
sim_10_filt <- runfixnames(100, sim_10_filt, correctnames)
jcp_unfilt <- runfixnames(100, jcp_unfilt, correctnames)
jcp_df <- runfixnames(100, jcp_df, correctnames)
jcp_sf <- runfixnames(100, jcp_sf, correctnames)
sg_uf <- runfixnames(100, sg_uf, correctnames)
sg_df <- runfixnames(100, sg_df, correctnames)
sg_sf <- runfixnames(100, sg_sf, correctnames)
sn_uf <- runfixnames(100, sn_uf, correctnames)
sn_df <- runfixnames(100, sn_df, correctnames)
sn_sf <- runfixnames(100, sn_sf, correctnames)
sim_unfilt_b$tip.label <- fixthenames(correctnames, sim_unfilt_b$tip.label)
sim_10_filt_b$tip.label <- fixthenames(correctnames, sim_10_filt_b$tip.label)
jcp_unfiltb$tip.label <- fixthenames(correctnames, jcp_unfiltb$tip.label)
jcp_dfb$tip.label <- fixthenames(correctnames, jcp_dfb$tip.label)
jcp_sfb$tip.label <- fixthenames(correctnames, jcp_sfb$tip.label)
sg_ufb$tip.label <- fixthenames(correctnames, sg_ufb$tip.label)
sg_dfb$tip.label <- fixthenames(correctnames, sg_dfb$tip.label)
sg_sfb$tip.label <- fixthenames(correctnames, sg_sfb$tip.label)
sn_ufb$tip.label <- fixthenames(correctnames, sn_ufb$tip.label)
sn_dfb$tip.label <- fixthenames(correctnames, sn_dfb$tip.label)
sn_sfb$tip.label <- fixthenames(correctnames, sn_sfb$tip.label)



all_trees <- c(sim_unfilt, sim_10_filt, 
               sim_unfilt_b, sim_10_filt_b,jcp_unfilt, jcp_df, jcp_sf, jcp_unfiltb,
               jcp_dfb, jcp_sfb, sg_uf, sg_df, sg_sf, sg_ufb, sg_dfb, sg_sfb,
               sn_uf, sn_df, sn_sf, sn_ufb, sn_dfb, sn_sfb)
class(all_trees) <- "multiPhylo"

names(all_trees)[1:100] <- paste0("sim_uf",1:100)
names(all_trees)[101:200] <- paste0("sim_df",1:100)
names(all_trees)[201]<- paste0("sim_ufb")
names(all_trees)[202]<- paste0("sim_dfb")
names(all_trees)[203:302] <- paste0("jcp_uf",1:100)
names(all_trees)[303:402] <- paste0("jcp_df",1:100)
names(all_trees)[403:502] <- paste0("jcp_sf",1:100)
names(all_trees)[503]<- paste0("jcp_ufb")
names(all_trees)[504]<- paste0("jcp_dfb")
names(all_trees)[505]<- paste0("jcp_sfb")
names(all_trees)[506:605] <- paste0("sg_uf",1:100)
names(all_trees)[606:705] <- paste0("sg_df",1:100)
names(all_trees)[706:805] <- paste0("sg_sf",1:100)
names(all_trees)[806]<- paste0("sg_ufb")
names(all_trees)[807]<- paste0("sg_dfb")
names(all_trees)[808]<- paste0("sg_sfb")
names(all_trees)[809:908] <- paste0("sn_uf",1:100)
names(all_trees)[909:1008] <- paste0("sn_df",1:100)
names(all_trees)[1009:1108] <- paste0("sn_sf",1:100)
names(all_trees)[1109]<- paste0("sn_ufb")
names(all_trees)[1110]<- paste0("sn_dfb")
names(all_trees)[1111]<- paste0("sn_sfb")

Dtype <- c(rep("Unfiltered",100),rep("Default",100),
           rep("Best Tree",2),rep("Unfiltered",100),rep("Default",100),
           rep("Standard",100),rep("Best Tree",3),rep("Unfiltered",100),rep("Default",100),
           rep("Standard",100),rep("Best Tree",3),rep("Unfiltered",100),rep("Default",100),
           rep("Standard",100),rep("Best Tree",3))

Dcols <- c("#999999","#56B4E9", "#009E73", "#E69F00")

Dscape <- treespace(all_trees,method='RF', nf=10)
symbs <- c(rep("Simulated",202),rep("JCP", 303),rep("SNiPgenie", 303),rep("snippy", 303)) 
plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Parameters",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             symbol_lab= "Pipeline",
             symbols= c("Simulated"="circle", "JCP"="triangle","SNiPgenie"="square","snippy"="star"),
             ellipses=TRUE)


fixthenames <- function(correct, query){
  
  outvector <- query
  
  for(position in correct){
   
    place <- grep(position, query)
    
    outvector[place] <- position
    
  }
  for(spot in 1:length(outvector)){
      if(outvector[spot]=="REF"|| outvector[spot] == "Reference" || outvector[spot] == "Ref-1997" || outvector[spot] == "ref"){
        
        outvector[spot] <- "Ref"
        
      }
      
      
  }
  return(outvector)
}

runfixnames <- function(num, into, correct){
  
  out <- into
  
  for(x in 1:num){
    
    out[[x]]$tip.label <- fixthenames(correct, into[[x]]$tip.label)
  }
  return(out)
}