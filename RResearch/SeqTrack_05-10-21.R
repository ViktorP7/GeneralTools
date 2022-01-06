# Running seqtrack on problem herds

# load packages
library("ape")
library("adegenet")
library("igraph")
library("ggplot2")
library("plyr")

# load the SNP sequence data
snp <-read.dna("C:/Users/UCD/Desktop/UbuntuSharedFolder/Cork10Fastqs/C10genie/nrc10core.fa", format="fa")
snpt <-read.dna("C:/Users/UCD/Desktop/UbuntuSharedFolder/SnipGenieTyr/nrtcore.fa", format="fa")



# load meta
cases <- read.csv("C:/Users/UCD/Desktop/UbuntuSharedFolder/Cork10Fastqs/C10genie/MAP_cases.csv")
tyrcases <- read.csv("C:/Users/UCD/Desktop/UbuntuSharedFolder/SnipGenieTyr/Tyr_cases.csv")


# define the columms with dates as GMT dates in POSIXct format.
birth <- as.POSIXct(cases$birthd, tz = "GMT") 
tbirth <- as.POSIXct(tyrcases$birthd, tz = "GMT") 


# create distance matrix
distmat <- dist.dna(snp, model="N", pairwise.deletion = TRUE, as.matrix = TRUE)
tyrdistmat <- dist.dna(snpt, model="N", pairwise.deletion = TRUE, as.matrix = TRUE)
write.csv(distmat, "DistMat_pairwDeletion.csv")


# define length of nucleotide sequence and mutation rate mu 
nbNucl <- ncol(as.matrix(snp))
tnbNucl <- ncol(as.matrix(snpt))
mu_MAP <- 2.91e-08*2/365    ## 0.25 substitutions per genome per year

Date <- birth
tdate <-tbirth

# run SeqTrack analysis 
res_Basic <- seqTrack(distmat, x.names=cases$label, x.dates=Date, mu=mu_MAP, haplo.le=nbNucl)
tres_Basic <- seqTrack(tyrdistmat, x.names=tyrcases$label, x.dates=tdate, mu=mu_MAP, haplo.le=tnbNucl)

write.csv(res_Basic, file = "C10_summary.csv")
write.csv(tres_Basic, file = "Tyr_summary.csv")


# calculate statistical support for inferred ancestries 
p_Basic <- get.likelihood(res_Basic, mu=mu_MAP, haplo.length=nbNucl)
tp_Basic <- get.likelihood(tres_Basic, mu=mu_MAP, haplo.length=tnbNucl)


# replace all ancestors with weight > "x" with NA. In the article we used x = 6. 
x <- 6
res_Basic$ances[res_Basic$weight > x] <- NA
tres_Basic$ances[tres_Basic$weight > x] <- NA


# create summary files: "ares" contains ancestries, "ap" contains p-values of statistical support for ancestries  
ares <- data.frame(res_Basic$ances)
ap <- data.frame(p_Basic)


# save time stamp-specific summary files, begin with "birth"
ares_birth <- ares
ap_birth <- ap

# plot pairwise genomic distance of all snp sequences 
hist <- hist(distmat,  col="lightgrey", nclass=250, xlim = c(0,50), ylim = c(0,1200), 
             main="Distribution of pairwise genomic distances (Cork 10)",
             xlab="Number of differing SNPs")
hist <- hist(tyrdistmat,  col="lightgrey", nclass=250, xlim = c(0,40), ylim = c(0,200), 
             main="Distribution of pairwise genomic distances (Tyrone CaA)",
             xlab="Number of differing SNPs")

# plot SeqTrack transmission trees 
outputFile <- paste("C10_ST_06-10-21.pdf", sep="")
pdf(outputFile, height=50, width=50)

res_network <- plot(res_Basic)
mtext(side=3, text="dark red: no/few mutations; light grey: many mutations", cex = 10)
res_network

dev.off()

#Create interactive plot for SeqTrack tree
tkplot(res_network, vertex.color="lightblue", layout=layout_with_lgl, 
       vertex.size=13, vertex.label.cex = 1.5,
       edge.label.cex = 1, canvas.width = 500, canvas.height = 500)

tres_network <- plot(tres_Basic)

tkplot(tres_network, vertex.color="lightblue", layout=layout_with_lgl, 
       vertex.size=13, vertex.label.cex = 1.5,
       edge.label.cex = 1.5, canvas.width = 500, canvas.height = 500)


# 3.) ANALYSIS OF RECONSTRUCTED MAP TRANSMISSION TREES

# create table with number of descending isolates and number of offspring
a <- ares_birth
ancesfreq_Basic <- as.data.frame(table(a$res_Basic))

cas <- data.frame(cow=cases$cow, MAP_id=cases$id)
cas$freq_Basic = ancesfreq_Basic$Freq[match(cases$id, ancesfreq_Basic$Var1)]  # gives number of descending isolates at isolate level

cas[is.na(cas)] <- 0                                                          # replaces all "NA" by 0

cow_freq_Basic <- count(cas, "cow", "freq_Basic")                             # gives number of descending isolates at cow level



cow_freqb <- data.frame(cow=cow_freq_Basic$cow, b_Basic=cow_freq_Basic$freq)
cas_birth <- cas


cow_freqs <- data.frame(cow=cow_freq_Basic$cow, s_Basic=cow_freq_Basic$freq)

# save results for data handling and further analysis in R or in another software
# data handling outside R: accout for potential right and left censoring of data by using 
# a subsample of all isolates to investigate the role of an individual cow in MAP spread (see 4.) 
cow_freqbs <- data.frame(cow_freqb, cow_freqs)
cas_bs <- data.frame(b=cas_birth)

write.csv(cas_bs, "cas_bs.csv")
write.csv(cow_freqbs, "cow_freqbs.csv")


# 4.) CORRELATION BETWEEN RECONSTRUCTED NUMBER OF OFFSPRING AND DISEASE PHENOTYPES

## load cow disease phenotypes and concatenate them to final (adapted) table with number of offspring per cow.
# format of table with cow disease phenotypes: "cow_id", "disease phenotypes" in binary or categorical coding.
# make sure that both tables are sorted in the same order of cows.
pheno <- read.csv("C:/Users/UCD/Desktop/UbuntuSharedFolder/Cork10Fastqs/C10genie/sheds.csv")
cow_freqbsx <- cow_freqbs
phenodf <- data.frame(pheno, cow_freqbsx)
write.csv(phenodf, "Cow_freqbsx_phenodf.csv")

## multifactorial analysis: correlate number of offspring (by scenario) with disease phenotype
# categorical coded phenotype - example: "shedding level"
cor.test(phenodf$b_Basic,phenodf$Sheddinglevel, method="pearson")


# boxplots with numbers of offspring produced by individual cows, by disease phenotype and scenario (Fig. 5 in article)
boxplot(phenodf$b_Basic~phenodf$Sheddinglevel,
        main="Shedding Level",
        ylab="Transmission Network Offspring (n)",
        xlab = "Shedding Level Cork 10")
legend("topleft", legend = c("0 - Unknown", "1 - High", "2 - Mid", "3 - Low"), title = "Shedding Levels")
