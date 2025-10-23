# Load required packages
library(venn)
library(VennDiagram)
library(ggplot2)
library(dplyr)

setwd("ADblood_lnc/ML/")
dir.results <- "vennresults"
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)


# Create feature list for Venn diagram
venn_list <- list(
  Xgboost = result_xgboost,
  RF = result_rfFuncs,
  Boruta = result_boruta,
  LassoLR = result_lasso,
  Pamr = result_pamr
)

# Basic Venn diagram
venn_basic <- venn(
  venn_list,
  ilabels = "counts",
  zcolor = 'style',  # Use default color style
  opacity = 0.3,     # Adjust color transparency
  box = FALSE,       # No border
  ilcs = 1.5,        # Number size
  sncs = 1.2         # Group name font size
)