# Load required packages
library(pROC)
library(ggplot2)
library(dplyr)
library(caret)

setwd("ADblood_lnc/ModelValidation/")
dir.results <- "results"
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

data_ADNI <- readRDS("/ADblood_lnc/datasets/data_ADNI.RDS")
data_AIBL <- readRDS("/ADblood_lnc/datasets/data_AIBL.RDS")
data_AddNeuroMed <- readRDS("/ADblood_lnc/datasets/data_AddNeuroMed.RDS")

preprocess_dataset <- function(data, dataset_name) {
  cat("Processing", dataset_name, "dataset...\n")
  
  # Factor conversion
  data$clin <- factor(data$clin, levels = c("CN", "AD"), labels = c("0", "1"))
  data$age <- as.numeric(data$age)
  
  # Data validation
  cat("  Dimensions:", dim(data), "\n")
  cat("  Clinical outcome distribution:\n")
  print(table(data$clin))
  cat("  Age range:", range(data$age, na.rm = TRUE), "\n")
  
  return(data)
}

# Preprocess all datasets
data_ADNI <- preprocess_dataset(data_ADNI, "ADNI")
data_AIBL <- preprocess_dataset(data_AIBL, "AIBL")
data_AddNeuroMed <- preprocess_dataset(data_AddNeuroMed, "AddNeuroMed")

# Define significant lncRNAs
significant_lncRNAs <- Reduce(intersect(result_boruta,result_pamr,result_xgboost,result_lasso,result_rfFuncs))

# Create formula for the model
model_formula <- as.formula(paste("clin ~", paste(significant_lncRNAs, collapse = " + ")))

model <- glm(
  formula = model_formula,
  data = data_AIBL,
  family = binomial(link = 'logit'),
  control = list(maxit = 100)
)

summary(model)
car::vif(model)



# Function for ROC analysis and visualization
perform_roc_analysis <- function(model, test_data, dataset_name, plot_title) {
  cat("Performing ROC analysis for", dataset_name, "...\n")
  
  # Generate predictions
  test_pred <- predict(model, newdata = test_data, type = "link")
  
  # Calculate ROC curve and AUC
  test_roc <- roc(test_data$clin, test_pred)
  auc_value <- auc(test_roc)
  
  cat("  AUC for", dataset_name, ":", round(auc_value, 4), "\n")
  
  # Create ROC plot
  roc_plot <- ggroc(test_roc, color = "#D60047FF", size = 1.2) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = plot_title,
      subtitle = paste("AUC =", round(auc_value, 4)),
      x = "Specificity",
      y = "Sensitivity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    annotate("text", x = 0.7, y = 0.3, 
             label = paste("AUC =", round(auc_value, 4)),
             size = 5, color = "black")
    print(roc_plot)
  return(list(roc = test_roc, auc = auc_value, plot = roc_plot))
}


roc_AIBL <- perform_roc_analysis(model, data_AIBL, "AIBL", "Training Set: AIBL (LncRNA Model)")
roc_ADNI <- perform_roc_analysis(model, data_ADNI, "ADNI", "Validation Set: ADNI (LncRNA Model)")
roc_AddNeuroMed <- perform_roc_analysis(model, data_AddNeuroMed, "AddNeuroMed", "Validation Set: AddNeuroMed (LncRNA Model)")


# Function to calculate AD scores
calculate_AD_scores <- function(data, model, dataset_name) {
  cat("Calculating AD scores for", dataset_name, "...\n")
  
  # Extract coefficients
  coefficients <- as.data.frame(coef(model))
  intercept_value <- coefficients[1, 1]
  
  # Select significant lncRNAs and ensure correct order
  significant_data <- data[, colnames(data) %in% significant_lncRNAs, drop = FALSE]
  significant_data <- significant_data[, match(rownames(coefficients)[-1], colnames(significant_data))]
  
  # Calculate scores
  scores <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    score <- intercept_value
    for (k in 1:ncol(significant_data)) {
      score <- score + significant_data[i, k] * coefficients[k + 1, 1]
    }
    scores[i] <- score
  }
  
  # Create results dataframe
  result <- data.frame(
    ID = rownames(data),
    score = scores,
    clin = data$clin,
    stringsAsFactors = FALSE
  )
  
  # Add original data
  result <- cbind(result, data)
  
  return(result)
}


score_AIBL <- calculate_AD_scores(data_AIBL, model, "AIBL")
score_ADNI <- calculate_AD_scores(data_ADNI, model, "ADNI")
score_AddNeuroMed <- calculate_AD_scores(data_AddNeuroMed, model, "AddNeuroMed")


# Save AD scores
saveRDS(score_AIBL, file.path(dir.results, "AD_scores_AIBL.RDS"))
saveRDS(score_ADNI, file.path(dir.results, "AD_scores_ADNI.RDS"))
saveRDS(score_AddNeuroMed, file.path(dir.results, "AD_scores_AddNeuroMed.RDS"))


# Train age + sex model on AIBL
model_age_sex <- glm(
  clin ~ age + sex,
  data = score_AIBL,
  family = binomial(link = 'logit'),
  control = list(maxit = 100)
)


summary(model_age_sex)


# ROC analysis for age + sex model
roc_age_sex_AIBL <- perform_roc_analysis(model_age_sex, score_AIBL, "AIBL", "AIBL: Age + Sex Model")
roc_age_sex_ADNI <- perform_roc_analysis(model_age_sex, score_ADNI, "ADNI", "ADNI: Age + Sex Model")
roc_age_sex_AddNeuroMed <- perform_roc_analysis(model_age_sex, score_AddNeuroMed, "AddNeuroMed", "AddNeuroMed: Age + Sex Model")


# Create performance comparison
performance_comparison <- data.frame(
  Dataset = rep(c("AIBL", "ADNI", "AddNeuroMed"), 2),
  Model = rep(c("LncRNA", "Age+Sex"), each = 3),
  AUC = c(
    roc_AIBL$auc, roc_ADNI$auc, roc_AddNeuroMed$auc,
    roc_age_sex_AIBL$auc, roc_age_sex_ADNI$auc, roc_age_sex_AddNeuroMed$auc
  )
)


# Create performance comparison plot
performance_plot <- ggplot(performance_comparison, aes(x = Dataset, y = AUC, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = round(AUC, 3)), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("LncRNA" = "#D60047FF", "Age+Sex" = "#00A08AFF")) +
  labs(
    title = "Model Performance Comparison Across Datasets",
    x = "Dataset",
    y = "AUC",
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  ylim(0, 1)

print(performance_plot)