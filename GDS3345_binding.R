## Getting data from GEO
## ==

library(GEOquery)
library(readr)

# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")

GDS3345 <- getGEO("GDS3345", GSEMatrix=FALSE, AnnotGPL=FALSE, getGPL=TRUE)
GDS3345_eset <- GDS2eSet(GDS3345)
GDS3345_eset_exprs <- exprs(GDS3345_eset)

## Renaming "Gene_Symbol" Column
## ==
colnames(fData(GDS3345_eset))[3] <- "Gene_Symbol"
GDS3345_Gene_Symbol <- fData(GDS3345_eset)[3]

## Saving as data frame
## ==
GDS3345_df <- data.frame(GDS3345_eset_exprs)

#Combining "Gene_Symbol" Column with Expression Data

GDS3345_bound <- cbind(GDS3345_Gene_Symbol, GDS3345_df) # 12625 obs

# Removing NAs
# ==

GDS3345_bound <- GDS3345_bound [!(is.na(GDS3345_bound$Gene_Symbol) | GDS3345_bound$Gene_Symbol == ""),]
  # 12181 obs

uniq_genesymb_GDS3345 <- unique(GDS3345_bound$Gene_Symbol)
data_list <- as.list(c())

for (i in uniq_genesymb_GDS3345)
{
  data_subset <- subset(GDS3345_bound, Gene_Symbol == i, select = -c(1:2))
  mat_subset <- as.matrix(data_subset)
  data_list[[i]] <- mat_subset
}

length(data_list) # [1] 9155

means_list <- lapply(data_list, colMeans)
mean_df <- data.frame(t(sapply(means_list, c)))

write_csv(mean_df, path = "GDS3345_mean_df.csv", append = FALSE, col_names = TRUE)