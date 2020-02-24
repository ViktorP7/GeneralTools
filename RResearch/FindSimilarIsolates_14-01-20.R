# Script to make summary of SNP distances between isolates
# Requires the EU tree to be generated in previous scripts

# Load packages
library(ape)

# Calculate genetic distance between isolates
euMat <- cophenetic(realEU)

# Round the distances
euMat <- round(euMat)

# Remove unecessary zeroes by filling with NA
for(index in 1:nrow(euMat)){
  
  euMat[index, index] <- NA
}

# Set the duplicate half of the matrix to NA
euMat[upper.tri(euMat)] <- NA

# Vector with all county names for the map
counties <- c("Antrim","Armagh","Carlow","Cavan","Clare","Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare","Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo","Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone","Waterford","Westmeath","Wexford","Wicklow")

# Run function to find similar isolates for group A
similarsA <- findSimilars(euroAmat, 30, counties)

# Run function to find similars for group F
similarsF <- findSimilars(euroFmat, 30, counties)

# Group G
similarsG <- findSimilars(euroGmat, 30, counties)

# Group C
similarsC <- findSimilars(euroCmat, 30, counties)

# Group B
similarsB <- findSimilars(euroBmat, 30, counties)

# Group D
similarsD <- findSimilars(euroDmat, 30, counties)

# Function to loop thru the matrix and find locations within the SNP threshold
findSimilars <- function(mat, threshold, counties){
  
  # Create result vector
  vector <- c()
  
  # Fetch the names for the matrix
  matNames <- rownames(mat)
  
  # Loop thru each row 
  for(row in 1:nrow(mat)){
    
    # Loop thru each col
    for(col in 1:ncol(mat)){
      
      # Skip yokey if it's NA
      if(is.na(mat[row,col])){
        
        next
      } else if(mat[row,col] > threshold){
        
        next
      } else{
        
        # Split to check location
        loc1 <- strsplit(matNames[row], split = "_")[[1]][2]
        loc2 <- strsplit(matNames[col], split = "_")[[1]][2]
        
        loc3 <- strsplit(matNames[row], split = "_")[[1]][1]
        loc4 <- strsplit(matNames[col], split = "_")[[1]][1]
        
        if(loc1 %in% counties ==  TRUE && loc2 %in% counties == TRUE){
          
          next
        } else if(loc3 == loc4){
          
          
          next
        }else if(loc1 %in% counties ==  TRUE && grepl("Ireland", loc4) == TRUE){
          
          
          next
        }else if(loc2 %in% counties ==  TRUE && grepl("Ireland", loc3) == TRUE){
          
          
          next
        } else { 
          # Paste together the row and col names and SNP count
          current <- paste(matNames[row], "~", matNames[col], "~", mat[row,col])
          
          vector <- append(vector, current)
        }
      }
    }
  }
  return(vector)
}
 