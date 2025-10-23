# Load required packages
library(rms)
library(pROC)
library(ggplot2)
library(dplyr)
library(regplot)

setwd("ADblood_lnc/Nomogram/")
dir.results <- "results"
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

score_AIBL <- readRDS("ADblood_lnc/ModelValidation/results/AD_scores_AIBL.RDS")
score_ADNI <- readRDS("ADblood_lnc/ModelValidation/results/AD_scores_ADNI.RDS")
score_AddNeuroMed <- readRDS("ADblood_lnc/ModelValidation/results/AD_scores_AddNeuroMed.RDS")

# Define model variables
model_vars <- c("clin", "score", "sex", "age", "Bas", "Bmem", "Bnv", 
                "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
# Create model dataset from AIBL (training set)
score_AIBL_model <- score_AIBL[, model_vars]

# Ensure immunity variables are numeric
immunity_vars <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", 
                   "Eos", "Mono", "Neu", "NK", "Treg")

score_AIBL_model[immunity_vars] <- lapply(score_AIBL_model[immunity_vars], as.numeric)

# Ensure categorical variables are factors
score_AIBL_model$sex <- factor(score_AIBL_model$sex)
score_AIBL_model$clin <- factor(score_AIBL_model$clin)
str(score_AIBL_model)

# Set up datadist for rms package
ddist <- datadist(score_AIBL_model)

# Manually set continuous variable ranges
for (var in immunity_vars) {
  ddist$limits[var, "Low:effect"] <- min(score_AIBL_model[[var]], na.rm = TRUE)
  ddist$limits[var, "High:effect"] <- max(score_AIBL_model[[var]], na.rm = TRUE)
  ddist$limits[var, "Adjust to"] <- median(score_AIBL_model[[var]], na.rm = TRUE)
}

# Set global datadist options
options(datadist = "ddist")

# Fit logistic regression model for nomogram
model_nomogram <- lrm(clin ~ score + sex + age + Bmem + Treg, 
                      data = score_AIBL_model, x = TRUE, y = TRUE)

nom <- nomogram(model_nomogram, 
                fun = plogis,
                fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                funlabel = "Predicted Probability",
                lp = FALSE,
                abbrev = FALSE)



CHDfit <- lrm(clin ~ score + Bmem + Treg, 
                              data = score_AIBL_model, x = TRUE, y = TRUE)
regplot(CHDfit, 
        plots=c("density","boxes") ,
        observation=FALSE,
        odds=FALSE,
        showP=FALSE,
        points=TRUE,
        clickable=FALSE,
        interval="confidence",
        dencol="#bebebe",
        boxcol="#bebebe")

