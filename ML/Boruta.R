##Boruta
# Load required packages
#install.packages("Boruta")
library(Boruta)
library(ggplot2)
library(dplyr)
setwd("/ADblood_lnc/ML/Boruta/")
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

# Perform Boruta feature selection
set.seed(20241218)

boruta.train <- Boruta(
  clin ~ .,
  data = data,
  doTrace = 2,
  maxRuns = 100,
  getImp = getImpRfZ  # Using random forest with z-score
)
print(boruta.train)

# Handle tentative features
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)

# Extract confirmed features
result_boruta <- as.data.frame(attStats(final.boruta))
confirmed_features <- rownames(result_boruta)[result_boruta$decision == "Confirmed"]
cat("\nNumber of confirmed features:", length(confirmed_features), "\n")
cat("Confirmed features:", paste(confirmed_features, collapse = ", "), "\n")

saveRDS(confirmed_features, file.path(dir.results, "confirmed_features_boruta.RDS"))

# Create visualization
boruta_plot <- plot(final.boruta, 
                    xlab = "Features", 
                    xaxt = "n",
                    main = "Boruta Feature Importance")

#ggplot2
feature_importance <- attStats(final.boruta)
feature_importance$feature <- rownames(feature_importance)
feature_importance$decision <- factor(feature_importance$decision,
                                     levels = c("Confirmed", "Tentative", "Rejected"))
ggplot_boruta <- ggplot(feature_importance, 
                       aes(x = reorder(feature, meanImp), 
                           y = meanImp, 
                           fill = decision)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Boruta Feature Importance",
    x = "Features",
    y = "Mean Importance",
    fill = "Decision"
  ) +
  scale_fill_manual(values = c("Confirmed" = "#0E986F", 
                              "Tentative" = "#796CAD", 
                              "Rejected" = "#D65813")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )
ggplot_boruta











set.seed(20241218)
boruta.train <- Boruta(clin~., data = data, doTrace = 2,maxRuns = 100)
print(boruta.train)
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
plot(final.boruta, xlab = "", xaxt = "n")
result_boruta <- as.data.frame(attStats(final.boruta))
result_boruta <- rownames(result_boruta)[which(result_boruta$decision=="Confirmed")]