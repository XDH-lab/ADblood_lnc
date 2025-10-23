# Load required packages
library(pamr)
library(ggplot2)
library(dplyr)


setwd("ADblood_lnc/ML/Pamr/")
dir.results <- "results"
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)
data <- readRDS("/ADblood_lnc/datasets/consistent_lnc_AIBL.RDS")
data$clin <- factor(
  data$clin,
  levels = c("CN", "AD"),
  labels = c("0", "1")
)
cat("Dataset dimensions:", dim(data), "\n")
print(table(data$clin))

# Prepare data for PAMR
x <- as.matrix(data[, -1])
x <- t(x)
y <- data.frame(group = data[, 1])
rownames(y) <- rownames(data)
y_vector <- y[, 1]

mydata <- list(
  x = x,
  y = factor(y_vector),
  geneid = rownames(x),
  genenames = paste("g", rownames(x), sep = "")
)
cat("PAMR data structure created:\n")
cat("  Features:", nrow(mydata$x), "\n")
cat("  Samples:", ncol(mydata$x), "\n")
cat("  Classes:", levels(mydata$y), "\n")

# Train PAMR model
set.seed(20241218)
mytrain <- pamr.train(mydata)

# Perform cross-validation
mycv <- pamr.cv(mytrain, mydata)
pamr.plotcv(mycv)

# Find optimal threshold based on cross-validation error
pamr.confusion(mycv,1.3)## Create confusion matrix at optimal threshold
pamr.plotcvprob(mytrain,mydata,1.3)## Plot class probabilities at optimal threshold
pamr_result <- pamr.listgenes(mytrain, mydata, threshold=1.3)# Extract significant genes at optimal threshold
significant_genes <- pamr_result[, 1]
cat("Number of significant genes:", length(significant_genes), "\n")
cat("Significant genes:", paste(significant_genes, collapse = ", "), "\n")

saveRDS(significant_genes, file.path(dir.results, "significant_genes_pamr.RDS"))

