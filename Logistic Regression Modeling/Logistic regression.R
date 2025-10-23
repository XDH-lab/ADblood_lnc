# Load required packages
library(dplyr)
library(foreach)
library(doParallel)

# Set up parallel processing
registerDoParallel(cores = parallel::detectCores() - 1)

# Function to perform logistic regression analysis for a dataset
perform_logistic_analysis <- function(beta_matrix, clinical_data, dataset_name) {
  cat("Performing logistic regression analysis for", dataset_name, "\n")
  
  results <- foreach(i = 1:ncol(beta_matrix), .combine = "rbind", .packages = "dplyr") %dopar% {
    if (i %% 100 == 0) {
      cat("Processing", dataset_name, "- CpG", i, "of", ncol(beta_matrix), "\n")
    }
    
    # Prepare data for current CpG site
    current_data <- data.frame(
      beta = beta_matrix[, i],
      clinical_data
    )
    
    # Factor conversion
    current_data$clin <- factor(
      current_data$clin, 
      levels = c("CN", "AD"), 
      labels = c("0", "1")
    )
    current_data$sex <- factor(
      current_data$sex,
      levels = c("Female", "Male"), 
      labels = c("0", "1")
    )
    current_data$age <- as.numeric(current_data$age)
    
    # Fit logistic regression model
    model <- glm(
      clin ~ beta + age + sex + Bas + Bmem + Bnv + CD4mem + CD4nv + 
            CD8mem + CD8nv + Eos + Mono + Neu + NK + Treg,
      data = current_data,
      family = binomial
    )
    
    # Extract results
    model_summary <- summary(model)$coefficients
    beta_coef <- model_summary["beta", "Estimate"]
    beta_pvalue <- model_summary["beta", "Pr(>|z|)"]
    
    # Determine direction based on coefficient
    if (beta_coef > 0) {
      direction <- "++"  # Up-regulated in AD
    } else if (beta_coef < 0) {
      direction <- "--"  # Down-regulated in AD
    } else {
      direction <- "NA"  # No change (theoretically rare)
    }
    
    data.frame(
      cpg = colnames(beta_matrix)[i],
      p_value = beta_pvalue,
      direction = direction,
      coefficient = beta_coef,
      stringsAsFactors = FALSE
    )
  }
  
  rownames(results) <- results$cpg
  return(results)
}

# Perform analysis for each dataset
## AIBL dataset
p_values_AIBL <- perform_logistic_analysis(
  beta_matrix = beta_AIBL,
  clinical_data = clin_AIBL,
  dataset_name = "AIBL"
)

## ADNI dataset  
p_values_ADNI <- perform_logistic_analysis(
  beta_matrix = beta_ADNI,
  clinical_data = clin_ADNI, 
  dataset_name = "ADNI"
)

## AddNeuroMed dataset
p_values_AddNeuroMed <- perform_logistic_analysis(
  beta_matrix = beta_AddNeuroMed,
  clinical_data = clin_AddNeuroMed,
  dataset_name = "AddNeuroMed"
)

# Stop parallel cluster
stopImplicitCluster()

# Identify significant lncRNAs in AIBL with consistent direction across datasets
identify_consistent_lncRNAs <- function(pvals_aibl, pvals_adni, pvals_anm, significance_threshold = 0.05) {
  
  # Find significant lncRNAs in AIBL
  significant_lnc <- rownames(pvals_aibl)[pvals_aibl$p_value <= significance_threshold]
  cat("Number of significant lncRNAs in AIBL:", length(significant_lnc), "\n")
  
  if (length(significant_lnc) == 0) {
    warning("No significant lncRNAs found in AIBL dataset")
    return(character(0))
  }
  
  # Check direction consistency across datasets
  direction_aibl <- pvals_aibl[significant_lnc, "direction"]
  direction_adni <- pvals_adni[significant_lnc, "direction"]
  direction_anm <- pvals_anm[significant_lnc, "direction"]
  
  consistent_direction <- (direction_aibl == direction_adni) & (direction_aibl == direction_anm)
  consistent_lnc <- significant_lnc[consistent_direction & !is.na(consistent_direction)]
  
  cat("Number of lncRNAs with consistent direction across all datasets:", length(consistent_lnc), "\n")
  
  return(consistent_lnc)
}

# Get consistent lncRNAs
consistent_lnc <- identify_consistent_lncRNAs(
  pvals_aibl = p_values_AIBL,
  pvals_adni = p_values_ADNI, 
  pvals_anm = p_values_AddNeuroMed
)

# Prepare and save consistent lncRNA data from AIBL
if (length(consistent_lnc) > 0) {
  consistent_lnc_AIBL <- beta_AIBL[, colnames(beta_AIBL) %in% consistent_lnc, drop = FALSE]
  consistent_lnc_AIBL <- data.frame(
    clin = clin_AIBL$clin,
    consistent_lnc_AIBL,
    stringsAsFactors = FALSE
  )
  
  # Save results
  saveRDS(consistent_lnc_AIBL, "/ADblood_lnc/datasets/consistent_lnc_AIBL.RDS")
  cat("Consistent lncRNA data saved to /ADblood_lnc/datasets/consistent_lnc_AIBL.RDS\n")
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  cat("Total CpGs analyzed in AIBL:", nrow(p_values_AIBL), "\n")
  cat("Significant CpGs in AIBL (p < 0.05):", sum(p_values_AIBL$p_value < 0.05), "\n")
  cat("CpGs with consistent direction across datasets:", length(consistent_lnc), "\n")
} else {
  warning("No consistent lncRNAs found across datasets")
}

# Save results
data_ADNI <- cbind(beta_ADNI, clin_ADNI)
data_AIBL <- cbind(beta_AIBL, clin_AIBL)
data_AddNeuroMed <- cbind(beta_AddNeuroMed, clin_AddNeuroMed)
saveRDS(data_ADNI, "/ADblood_lnc/datasets/data_ADNI.RDS")
saveRDS(data_AIBL, "/ADblood_lnc/datasets/data_AIBL.RDS")
saveRDS(data_AddNeuroMed, "/ADblood_lnc/datasets/data_AddNeuroMed.RDS")



saveRDS(p_values_AIBL, "/ADblood_lnc/datasets/p_values_AIBL.RDS")
saveRDS(p_values_ADNI, "/ADblood_lnc/datasets/p_values_ADNI.RDS") 
saveRDS(p_values_AddNeuroMed, "/ADblood_lnc/datasets/p_values_AddNeuroMed.RDS")

