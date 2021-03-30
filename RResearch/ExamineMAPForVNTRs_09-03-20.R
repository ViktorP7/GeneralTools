# Script to determine VNTR sequences from MAP K10
# Created on 09/03/20

# Load packages
library(seqinr)

# Designate path to fasta file
pathFasta <- "C:/Users/UCD/Desktop/UbuntuSharedFolder/Winter2018MAPSequencing/MAP-FASTQs/vcfFiles/Bryantandus/Reference/MAPK10.fasta"

# Read in the fasta 
K10Fasta <- read.fasta(pathFasta)

sites <- list(

  # Check VNTR site 292
  site292 <- K10Fasta$NC_002944.2[3253590:3253889],

  # Check site X3
  siteX3 <- K10Fasta$NC_002944.2[4441875:4442070],

  # Check site 25
  site25 <- K10Fasta$NC_002944.2[3665598:3665947],

  # Check site 47
  site47 <- K10Fasta$NC_002944.2[4128604:4128821],

  # Check site 3
  site3 <- K10Fasta$NC_002944.2[131320:131527],

  # Check site 7
  site7 <- K10Fasta$NC_002944.2[3711417:3711619],

  # Check site 10
  site10 <- K10Fasta$NC_002944.2[4279553:4279855],

  # Check site 32 
  site32 <- K10Fasta$NC_002944.2[1125707:1126004]
)

# Make vector of breaks
breaks <- sort(c(3253590,3253889,4441875,4442070,3665598,3665947,4128604,4128821,131320,131527,3711417,3711619,4279553,4279855,1125707,1126004))

# Define repeats
oldrepeats <- list(
  r292 <- strsplit(tolower("GTCATCTGCGCCGCTCCTCCTCATCGCTGCGCTCTGCATCGTCGTCGGCGCGA"),"")[[1]],
  rX3 <- strsplit(tolower("CGCGCGTATGACGATGCAGAGCGTAGCGATGAGGAGGAATCGCGCCGATGACC"),"")[[1]],
  r25 <- strsplit(tolower("AATGCTTCGTCGCCGGGCTCCACCCCAATCACCCACTCCTGCGCATCCCGCTGCGCGG"),"")[[1]],
  r47 <- strsplit(tolower("CCCGCGCCCGGCGAGACGGGCGCCCTGATGAGA"),"")[[1]],
  r3 <- strsplit(tolower("TGGCCGCGATGCCGACGCTGCCGGAGC"),"")[[1]],
  r7 <- strsplit(tolower("CGAAATATTCGCCGTGAGAACA"),"")[[1]],
  r10 <- strsplit(tolower("GATGGCGGCGGCCCGCCGCGCCCGGCTTCGCCGCGCTTGCGACCGCCGCTAGACC"),"")[[1]],
  r32 <- strsplit(tolower("GCCCGGGGGCGGGCCGTA"),"")[[1]]
)

repeats <- list(
  r292 <- strsplit(tolower("GTCATCTGCGCCGCTCCTCCTCATCGCTGCGCTCTGCATCGTCGTCGGCGCGA"),"")[[1]],
  rX3 <- strsplit(tolower("CGCGCGTATGACGATGCAGAGCGTAGCGATGAGGAGGAATCGCGCCGATGACC"),"")[[1]],
  r25 <- strsplit(tolower("CCCGCTGCGCGGAATGCTTCGTCGCCGGGCTCCACCCCAATC"),"")[[1]],
  r47 <- strsplit(tolower("GAGACGGGCGCCCTGATGAGACCCGCGCCCGG"),"")[[1]],
  r3 <- strsplit(tolower("GCTGGCCGCGATGCCGACGCTGCCGGA"),"")[[1]],
  r7 <- strsplit(tolower("CGAAATATTCGCCGTGAGAACA"),"")[[1]],
  r10 <- strsplit(tolower("GATGGCGGCGGCCCGCCGCGCCCGGCTTCGCCGCGCTTGCGACCGCCGCTAG"),"")[[1]],
  r32 <- strsplit(tolower("GGCGTGCCGTAGCCCGG"),"")[[1]]
)

cn1 <- makeRepeatGenome(sites, breaks, repeats, K10Fasta$NC_002944.2, 1)
cn2 <- makeRepeatGenome(sites, breaks, repeats, K10Fasta$NC_002944.2, 2)
cn3 <- makeRepeatGenome(sites, breaks, repeats, K10Fasta$NC_002944.2, 3)
cn4 <- makeRepeatGenome(sites, breaks, repeats, K10Fasta$NC_002944.2, 4)
cn5 <- makeRepeatGenome(sites, breaks, repeats, K10Fasta$NC_002944.2, 5)

# Define function to create genome with x amount of vntr repeats at each locus
makeRepeatGenome <- function(sites, breaks, repeats, genome, x){
  
  # Make segments of genome
  segment <- list(
    pre3 <- genome[1:breaks[1]],
    pre32 <- genome[breaks[2]:breaks[3]],
    pre292 <- genome[breaks[4]:breaks[5]],
    pre25 <- genome[breaks[6]:breaks[7]],
    pre7 <- genome[breaks[8]:breaks[9]],
    pre47 <- genome[breaks[10]:breaks[11]],
    pre10 <- genome[breaks[12]:breaks[13]],
    prex3 <- genome[breaks[14]:breaks[15]],
    postx3 <- genome[breaks[16]:length(genome)]
  )
  
  # Paste together all segments
  for(index in 1:length(segment)){
    segment[[index]] <- paste(segment[[index]], sep="", collapse = "")
  }
  # Create a list to store repeated regions
  regions <- list()
  
  # Loop thru the sites and create multiple copy sites
  for(index in 1:length(sites)){
    
    torepeat <- rep(repeats[[index]], x)
    
    start <- append(sites[[index]][1:20], torepeat)
    ending <- append(start, sites[[index]][length(sites[[index]])-20:length(sites[[index]])])
    
    regions[[index]] <- ending
  }
  
  # Collapse regions
  for(index in 1:length(regions)){
    regions[[index]] <- paste(regions[[index]], sep = "", collapse = "")
  }
  
  # Create the final genome
  finalGenome <- strsplit(trimws(paste(segment[[1]],regions[[5]],segment[[2]], regions[[8]], segment[[3]], regions[[1]], segment[[4]], 
                       regions[[3]], segment[[5]], regions[[6]], segment[[6]], regions[[4]], segment[[7]], regions[[7]],
                       segment[[8]], regions[[2]], segment[[9]])),"")[[1]]
  
  finalGenome <- finalGenome[finalGenome != " "]
  
  write.fasta(finalGenome, paste("K10","copynum",x, sep = "-"), file.out = paste("K10","copynum",x,".fasta", sep = "-"))
  
  return(finalGenome)
}

