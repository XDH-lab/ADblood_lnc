# Load required packages
library(caret)
library(randomForest)
library(ggplot2)
library(dplyr)

setwd("/ADblood_lnc/ML/RandomForest/")
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

# Perform recursive feature elimination (RFE) with Random Forest
set.seed(20241218)

# Define RFE control parameters
control <- rfeControl(
  functions = rfFuncs,           # Use Random Forest
  method = "LGOCV",              # Leave-group-out cross-validation
  number = 10,                   # 10-fold cross-validation
  verbose = TRUE,                # Show progress
  returnResamp = "final",        # Return final resampling results
  saveDetails = TRUE             # Save detailed results
)

# Perform RFE
candidate_genes <- colnames(data)
results <- rfe(
  x = data[, -1],                # Features (all columns except first)
  y = data$clin,                 # Clinical outcome
  metric = "Accuracy",           # Optimization metric
  sizes = 1:(length(candidate_genes) - 2), # Step size of 1
  rfeControl = control
)


# Extract final selected genes
final_genes <- predictors(results)
cat("Number of final selected genes:", length(final_genes), "\n")
cat("Final selected genes:", paste(final_genes, collapse = ", "), "\n")
saveRDS(final_genes, file.path(dir.results, "final_selected_genes_rf.RDS"))


accuracy_results <- results$results
optimal_index <- which.max(accuracy_results$Accuracy)
optimal_features <- accuracy_results$Variables[optimal_index]
optimal_accuracy <- accuracy_results$Accuracy[optimal_index]

cat("Optimal number of features:", optimal_features, "\n")
cat("Optimal accuracy:", optimal_accuracy, "\n")


# Accuracy vs Number of Features plot using base R (matching original style)
pdf(file = file.path(dir.results, "randomforest_accuracy_plot.pdf"), 
    width = 6, height = 4.5)

# Set plotting parameters
par(bty = "o", 
    mgp = c(2, 0.5, 0), 
    mar = c(3.1, 4.1, 2.1, 2.1),
    tcl = -0.25, 
    las = 1)

# Create the accuracy plot
plot(accuracy_results$Variables,
     accuracy_results$Accuracy,
     ylab = "",
     xlab = "Number of Features",
     col = "steelblue",
     pch = 16,
     main = "Random Forest Feature Selection")

# Add connecting lines
lines(accuracy_results$Variables, accuracy_results$Accuracy, col = "steelblue")

# Highlight optimal point
points(optimal_features, optimal_accuracy,
       col = "red",
       pch = 19,
       cex = 1.5)

# Add Y-axis label
mtext("Accuracy (Repeated Cross-Validation)", side = 2, line = 2.5, las = 3)

# Add arrow and text annotation (using base R arrows instead of shape::Arrows)
arrows(x0 = optimal_features - 7, 
       x1 = optimal_features - 2,
       y0 = optimal_accuracy, 
       y1 = optimal_accuracy,
       length = 0.1,    # Arrow head length
       angle = 20,      # Arrow head angle
       code = 2,        # Arrow head at end
       lwd = 2,
       col = "black")

# Add feature count information
text(x = optimal_features - 7,
     y = optimal_accuracy,
     labels = paste0("N = ", optimal_features),
     pos = 2,
     cex = 0.9)

dev.off()