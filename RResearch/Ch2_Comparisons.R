library(ape)
library(phytools)
library(scales)
library(VennDiagram)
library("nVennR")
library("qpcR")
library("grImport2")
library(rsvg)
library(sciplot)
library(Biostrings)


# Read in tables
sim_uf_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sim_uf.tsv",
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors=FALSE, 
                           check.names=FALSE)
sim_df_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sim_dfupd.tsv",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sg_uf_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sg_uf.txt",
                          header = TRUE,
                          sep = " ",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sg_df_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sg_df.txt",
                          header = TRUE,
                          sep = " ",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sg_sf_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sg_sf.txt",
                          header = TRUE,
                          sep = " ",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sn_uf_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sn_uf.tab",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sn_df_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sn_df.tab",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)
sn_sf_snps <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/sn_sf.tab",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors=FALSE, 
                          check.names=FALSE)

# Read JCP files (different format, need merging)
jcp_uf_fa <- readDNAStringSet(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-seq-uf.fasta",
                          format="fasta", seek.first.rec = TRUE)
jcp_uf_pos <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-pos-uf.txt",
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors=FALSE, 
                         check.names=FALSE)
jcp_df_fa <- readDNAStringSet(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-seq-df.fasta",
                              format="fasta", seek.first.rec = TRUE)
jcp_df_pos <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-pos-df.txt",
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors=FALSE, 
                         check.names=FALSE)
jcp_sf_fa <- readDNAStringSet(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-seq_sf.fasta",
                              format="fasta", seek.first.rec = TRUE)
jcp_sf_pos <- read.table(file = "C:/Users/UCD/Desktop/UbuntuSharedFolder/PipelineSimulation/SnpTabs/jcp-pos_sf.txt",
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors=FALSE, 
                         check.names=FALSE)

# Convert JCP into standard format
jcp_uf_snps <- data.frame(matrix(nrow = nrow(jcp_uf_pos), ncol = length(jcp_uf_fa)))
jcp_df_snps <- data.frame(matrix(nrow = nrow(jcp_df_pos), ncol = length(jcp_df_fa)))
jcp_sf_snps <- data.frame(matrix(nrow = nrow(jcp_sf_pos), ncol = length(jcp_sf_fa)))

jcp_uf_snps[,1] <- jcp_uf_pos
jcp_df_snps[,1] <- jcp_df_pos
jcp_sf_snps[,1] <- jcp_sf_pos

for(col in 2:ncol(jcp_uf_snps)){
  
  jcp_uf_snps[,col] <- strsplit(as.character(jcp_uf_fa), split= "")[[col-1]]
}

for(col in 2:ncol(jcp_df_snps)){
  
  jcp_df_snps[,col] <- strsplit(as.character(jcp_df_fa), split= "")[[col-1]]
}

for(col in 2:ncol(jcp_sf_snps)){
  
  jcp_sf_snps[,col] <- strsplit(as.character(jcp_sf_fa), split= "")[[col-1]]
}

# Remove unnecessary cols where needed
sg_uf_snps <- sg_uf_snps[,-2]
sg_df_snps <- sg_df_snps[,-2]
sg_sf_snps <- sg_sf_snps[,-2]
sn_uf_snps <- sn_uf_snps[,-c(1,3)]
sn_df_snps <- sn_df_snps[,-c(1,3)]
sn_sf_snps <- sn_sf_snps[,-c(1,3)]


sim_vs_sg_uf <- calculateSNPstats(sim_uf_snps, sg_uf_snps)
sim_vs_sg_df <- calculateSNPstats(sim_df_snps, sg_df_snps)
sim_vs_sg_sf <- calculateSNPstats(sim_df_snps, sg_sf_snps)
sim_vs_sn_uf <- calculateSNPstats(sim_uf_snps, sn_uf_snps)
sim_vs_sn_df <- calculateSNPstats(sim_df_snps, sn_df_snps)
sim_vs_sn_sf <- calculateSNPstats(sim_df_snps, sn_sf_snps)
sim_vs_jcp_uf <- calculateSNPstats(sim_uf_snps, jcp_uf_snps)
sim_vs_jcp_df <- calculateSNPstats(sim_df_snps, jcp_df_snps)
sim_vs_jcp_sf <- calculateSNPstats(sim_df_snps, jcp_sf_snps)

# Plot out specificities and precisions
boxplot(sim_vs_sg_uf$Specificity,
        #sim_vs_sg_df$Specificity,
        #sim_vs_sg_sf$Specificity,
        sim_vs_sn_uf$Specificity,
        #sim_vs_sn_df$Specificity,
        #sim_vs_sn_sf$Specificity,
        sim_vs_jcp_uf$Specificity,
        #sim_vs_jcp_df$Specificity,
        #sim_vs_jcp_sf$Specificity, 
        main = "Specificity of Pipelines", 
        names = c("SNiPGenie", "Snippy", "JCP"),
        ylab = "Specificity", las=1)
stripchart(sim_vs_sg_uf$Specificity, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
           pch = 4)
#stripchart(sim_vs_sg_df$Specificity, add = TRUE, at =2, 
#           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
#           pch = 4)
#stripchart(sim_vs_sg_sf$Specificity, add = TRUE, at =3, 
#           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
#           pch = 4)
stripchart(sim_vs_sn_uf$Specificity, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
           pch = 4)
#stripchart(sim_vs_sn_df$Specificity, add = TRUE, at =5, 
#           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
#           pch = 4)
#stripchart(sim_vs_sn_sf$Specificity, add = TRUE, at =6, 
#           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
#           pch = 4)
stripchart(sim_vs_jcp_uf$Specificity, add = TRUE, at =3, 
           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
           pch = 4)
#stripchart(sim_vs_jcp_df$Specificity, add = TRUE, at =8, 
#           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
#           pch = 4)
#stripchart(sim_vs_jcp_sf$Specificity, add = TRUE, at =9, 
#           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
#           pch = 4)
#legend("topleft", legend = c("SG-SNiPgenie", "SN-snippy", "JCP", "UF- Unfiltered", "DF-Default", "SF-Standard"), 
#       text.col = c("blue","red","green","black","black","black"),
#       bty = "n")

# Plot out specificities and precisions
boxplot(sim_vs_sg_uf$Precision,
        #sim_vs_sg_df$Precision,
        #sim_vs_sg_sf$Precision,
        sim_vs_sn_uf$Precision,
        #sim_vs_sn_df$Precision,
        #sim_vs_sn_sf$Precision,
        sim_vs_jcp_uf$Precision,
        #sim_vs_jcp_df$Precision,
        #sim_vs_jcp_sf$Precision, 
        main = "Precision of Pipelines", 
        names = c("SNiPGenie", "Snippy", "JCP"),
        ylab = "Precision", las=1, ylim=c(0.9,1))
stripchart(sim_vs_sg_uf$Precision, add = TRUE, at =1, 
           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
           pch = 4)
#stripchart(sim_vs_sg_df$Precision, add = TRUE, at =2, 
#           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
#           pch = 4)
#stripchart(sim_vs_sg_sf$Precision, add = TRUE, at =3, 
#           method = "jitter", vertical = TRUE, col = alpha("blue",0.5),
#           pch = 4)
stripchart(sim_vs_sn_uf$Precision, add = TRUE, at =2, 
           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
           pch = 4)
#stripchart(sim_vs_sn_df$Precision, add = TRUE, at =5, 
#           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
#           pch = 4)
#stripchart(sim_vs_sn_sf$Precision, add = TRUE, at =6, 
#           method = "jitter", vertical = TRUE, col = alpha("red",0.5),
#           pch = 4)
stripchart(sim_vs_jcp_uf$Precision, add = TRUE, at =3, 
           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
           pch = 4)
#stripchart(sim_vs_jcp_df$Precision, add = TRUE, at =8, 
#           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
#           pch = 4)
#stripchart(sim_vs_jcp_sf$Precision, add = TRUE, at =9, 
#           method = "jitter", vertical = TRUE, col = alpha("green",0.5),
#           pch = 4)
#legend("topright", legend = c("SG-SNiPgenie", "SN-snippy", "JCP", "UF- Unfiltered", "DF-Default", "SF-Standard"), 
#       text.col = c("blue","red","green","black","black","black"),
#       bty = "n")

mean(sim_vs_sg_uf$Specificity)
mean(sim_vs_sg_df$Specificity)
mean(sim_vs_sg_sf$Specificity)
mean(sim_vs_sn_uf$Specificity)
mean(sim_vs_sn_df$Specificity)
mean(sim_vs_sn_sf$Specificity)
mean(sim_vs_jcp_uf$Specificity)
mean(sim_vs_jcp_df$Specificity)
mean(sim_vs_jcp_sf$Specificity)

mean(sim_vs_sg_uf$Precision)
mean(sim_vs_sg_df$Precision)
mean(sim_vs_sg_sf$Precision)
mean(sim_vs_sn_uf$Precision)
mean(sim_vs_sn_df$Precision)
mean(sim_vs_sn_sf$Precision)
mean(sim_vs_jcp_uf$Precision)
mean(sim_vs_jcp_df$Precision)
mean(sim_vs_jcp_sf$Precision)
# Function to calculate precision and sensitivity
calculateSNPstats <- function(df1,df2){
  
  dfout <- data.frame(matrix(nrow = ncol(df1)-1, ncol = 2))
  colnames(dfout) <- c("Specificity", "Precision")
  
  toskip <- c(setdiff(df1[,1], df2[,1]), setdiff(df2[,1], df1[,1]))
  
  fane <- length(setdiff(df1[,1], df2[,1]))
  fapo <- length(setdiff(df2[,1], df1[,1]))
  
  df1 <- df1[-which(df1[,1] %in% toskip),]
  
  if(fapo != 0){
    df2 <- df2[-which(df2[,1] %in% toskip),]
  }
  for(col in 2:ncol(df1)){
     
    trpo <- 0
    
    for(row in 1:nrow(df1)){
        
      if(df1[row,col] == df2[row,col] || df2[row,col] == "N"){
          
        trpo <- trpo + 1
      } else {
          
        fane <- fane +1
      }
      
    }
    spec <- trpo/(trpo + fane)
    prec <- trpo/(trpo + fapo)
    
    dfout[col-1,1] <- spec
    dfout[col-1,2] <- prec
  }

  return(dfout)

}
