# Load required packages
library(xgboost)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

setwd("/ADblood_lnc/ML/XGBoost/")
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

# Prepare data for Xgboost
set.seed(20241218)
train_indices <- createDataPartition(
  data$clin,
  p = 0.7,
  list = FALSE
)

train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

train_x <- data.matrix(train_data[, -1])

test_x <- data.matrix(test_data[, -1])

# Create XGBoost DMatrix objects
xgb_train <- xgb.DMatrix(data = train_x, label = train_y)
xgb_test <- xgb.DMatrix(data = test_x, label = test_y)

# Define watchlist for monitoring training progress
watchlist <- list(train = xgb_train, test = xgb_test)

# Train initial model to determine optimal number of rounds
model_xgb_cv <- xgb.train(
  data = xgb_train,
  max.depth = 3,
  watchlist = watchlist,
  nrounds = 200
)

## Extract evaluation log to find optimal rounds ##

model_xgb_final <- xgboost(
  data = xgb_train,
  max.depth = 3,
  nrounds = optimal_rounds,##optimal_rounds
  verbose = 0
)

# Extract feature importance
feature_importance <- xgb.importance(model = model_xgb_final)
xgb.ggplot.importance(feature_importance)


result_xgboost <- feature_importance$Feature
saveRDS(result_xgboost, file.path(dir.results, "significant_genes_xgboost.RDS"))