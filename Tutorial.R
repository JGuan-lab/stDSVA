############################################################
# Tutorial: Run spatial transcriptomics deconvolution with stDSVA
# ------------------------------------------------------------
# This example demonstrates how to perform spatial deconvolution 
# using the stDSVA framework on a provided example dataset.
############################################################

# 1. Load the main stDSVA functions
source("stDSVA.R")  

# 2. Load the example dataset
# The dataset contains both single-cell reference data and ST data
example_dataset <- readRDS("example_dataset.rds")

# 3. Extract relevant matrices
# - scData: single-cell expression matrix (genes × cells)
# - meta: cell metadata, including cell type annotations
# - ST: spatial expression matrix (genes × spots)
scData <- as.matrix(example_dataset$comData$Train$data)
meta   <- data.frame(example_dataset$comData$Train$pData)
ST     <- as.matrix(example_dataset$spatial.real)

# 4. Initialize model parameters
# The getParas() function prepares parameter settings for model training.
# For this quick-start tutorial, we only use a single set of hyperparameters.
# In practice, you can specify vectors (e.g., lambdaC = c(0.1, 1, 100)) 
# to perform grid search tuning.
paras <- getParas(
  scData,
  meta,
  ST,
  windowTrain = 375,
  constraints_l = 1000,
  lambda1 = 0,
  lambda2 = 0,
  lambdaC = 1000
)

# 5. Run the deconvolution
# The runDeconv_my() function performs spatial deconvolution 
# using the parameter configuration obtained above.
res <- runDeconv_my(
  paras,
  constraints_l = 1000,
  lambda1 = 0,
  lambda2 = 0,
  lambdaC = 1000
)

# 6. Evaluate deconvolution performance
getPearsonRMSE(res$p, example_dataset$P_test)
getPearsonRMSE(res$c, example_dataset$C)

############################################################
# Notes:
# - getParas(): prepares parameter configurations for tuning.
#   For real applications, use multiple values for constraints_l/lambda1/lambda2/lambdaC
#   to perform grid search and select the best configuration.
#
# - runDeconv_my(): executes the deconvolution step using selected parameters.
#
############################################################
