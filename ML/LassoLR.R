# Load required packages
library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)

setwd("/ADblood_lnc/ML/LASSO/")
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


# Prepare data for Lasso
states <- as.matrix(data)
x <- states[, -1]  # Features (all columns except the first)
y <- as.numeric(states[, 1])  # Response variable (first column)

# Perform cross-validation to find optimal lambda
set.seed(20241218)

cvfit <- cv.glmnet(
  x = x,
  y = y,
  type.measure = "mse",
  nfolds = 10,
  alpha = 1,  # 1 for LASSO, 0 for Ridge
  family = "binomial"
)
# Extract optimal lambda values
lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
cat("Optimal lambda values:\n")
cat("  Lambda min:", lambda_min, "\n")
cat("  Lambda 1se:", lambda_1se, "\n")

# Create cross-validation plot
plot(cvfit, main = "LASSO Cross-Validation")

# Fit final LASSO model
lasso_model <- glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  nlambda = 100
)
print(lasso_model)

# L1 norm plot
plot(lasso_model, label = TRUE, main = "LASSO Coefficients vs L1 Norm")
# Lambda plot
plot(lasso_model, xvar = "lambda", label = TRUE, 
     main = "LASSO Coefficients vs Lambda")
# Deviance explained plot
plot(lasso_model, xvar = "dev", label = TRUE,
     main = "LASSO Coefficients vs Fraction Deviance Explained")

# Get coefficients for both lambda.min and lambda.1se
coef_matrix <- as.matrix(coef(lasso_model, s = c(lambda_min, lambda_1se)))
coef_df <- as.data.frame(coef_matrix)
colnames(coef_df) <- c("lambda_min", "lambda_1se")

# Identify non-zero coefficients
non_zero_min <- rownames(coef_df)[which(coef_df$lambda_min != 0)]
non_zero_1se <- rownames(coef_df)[which(coef_df$lambda_1se != 0)]



result_lasso_min <- setdiff(non_zero_min, "(Intercept)")
result_lasso_1se <- setdiff(non_zero_1se, "(Intercept)")

saveRDS(result_lasso_min, file.path(dir.results, "significant_genes_lasso_min.RDS"))