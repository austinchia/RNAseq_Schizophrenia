# 1) Split dataframe into 2 classes - disease and control
  # take from metadata

library(readr)
library("genefilter")

data <- read.csv("GDS3345_mean_df.csv", header = TRUE, stringsAsFactors =  FALSE, sep = ",")
rownames(data) <- rownames(mean_df)
data <- cbind(data[,1:12], data[,38:49])

# control class - 1-12
# disease class - 13-24

# 2) Genefilter package method
#BiocManager::install("genefilter")

library(progress)
pb <- progress_bar$new(total = 1000)
fac <- c(rep("normal",4),rep("disease",4))
data_mat <- as.matrix(data)
boot_list <- list()

for (i in 1:1000) {
  my_significant_genes <- c()
  control <- as.matrix(sample(data[,1:12], size = 4, replace = TRUE))
  disease <- as.matrix(sample(data[,13:24], size = 4, replace =TRUE))
  m3 <- cbind(control, disease)
  test_ttest <- rowttests(m3,fac = as.factor(fac))
  my_significant_genes_v2 <- as.numeric(test_ttest$p.value < 0.05)
  boot_list <- append(boot_list, list(my_significant_genes_v2))
  
  pb$tick()
  Sys.sleep(1 / 1000)
}

## 3) Transforming nested list into dataframe

boot_df <- data.frame(matrix(unlist(boot_list), nrow=nrow(data), byrow=F),stringsAsFactors=FALSE)
write_csv(boot_df, path = "GDS3345_boot_df.csv", append = FALSE, col_names = TRUE)

## 4) Jaccard Coefficient ##

jac_func <- function (x, y) {
  M_11 <- sum(x == 1 , y == 1)
  M_10 <- sum(x == 1 , y == 0)
  M_01 <- sum(x == 0 , y == 1)
  return (M_11 / (M_11 + M_10 + M_01))
}

jac_df <- data.frame(matrix(data = NA, nrow = length(boot_df), ncol = length(boot_df)))
pb <- progress_bar$new(total = 1000)

for (r in 1:length(boot_df)) {
  for (c in 1:length(boot_df)) {
    if (c == r) {
      jac_df[r,c] = 1
    } else if (c > r) {
      jac_df[r,c] = jac_func(boot_df[,r], boot_df[,c])
    }
  }
  pb$tick()
  Sys.sleep(1 / 1000)
}

variable_names <- sapply(boot_df, attr, "label")
colnames(jac_df) <- paste0("S", seq(1:ncol(jac_df)))
rownames(jac_df) <- paste0("S", seq(1:nrow(jac_df)))

jac_mat <- matrix(jac_df)


## Determining significance of genes

sig_hist <- hist(row_total)

# For data
row_total <- as.matrix(rowSums(boot_df))
sig_gene <- as.numeric(row_total>200)
dat_fac <- as.factor(sig_gene)

# For reference
ref_fac <- as.factor(boot_df[,1])

library(caret)
cm <- confusionMatrix(data = dat_fac, reference = ref_fac )

draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Class1', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Class2', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Class1', cex=1.2, srt=90)
  text(140, 335, 'Class2', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

dat_cm <- draw_confusion_matrix(cm)
