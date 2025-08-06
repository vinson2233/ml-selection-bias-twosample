# ===================================================================
# Synthetic Data Generator for Testing Multiple Datasets
# ===================================================================
# This script generates synthetic datasets for testing the 
# statistical estimation framework with different data structures.

library(dplyr)
library(tidyr)

# Set random seed for reproducibility
set.seed(123)

# ===================================================================
# Helper
# ===================================================================
add_noise <- function(signal,target_r2){
  noise_sd = sqrt(var(signal)*(1 - target_r2)/target_r2)
  y = signal + rnorm(length(signal),0,sd = noise_sd)
  return(y)
}

# ===================================================================
# Synthetic Dataset Generation
# ===================================================================
# Follow the procedure described by "Chen, Y., Li, P., & Wu, C. (2019). Doubly Robust Inference With Nonprobability Survey Samples. Journal of the American Statistical Association, 115(532), 2011-2021. https://doi.org/10.1080/01621459.2019.1677241"
#' Generate a simple synthetic dataset with known relationships
#' 
#' @param N Number of rows in the dataset
#' @param sigma The amount of noise to be introduced. in the true regression model, $\sigma$ chosen such that correlation between y and the linear prediction is 0.3, 0.5 and 0.8 (From manual testing,$\sigma$ = 10 gives cor 0.3. $\sigma$ = 5.5 gives cor 0.5,$\sigma$ = 2.4 gives cor 0.8)  

generate_synthetic_linear_low_noise <- function(N = 1000000) {
  z1 <- rbinom(N, 1, 0.5)
  z2 <- runif(N, 0, 2)
  z3 <- rexp(N, rate=1)
  z4 <- rchisq(N, df=4)
  
  x1 = z1
  x2 = z2 + 0.3*x1
  x3 = z3 + 0.2*(x1 + x2)
  x4 = z4 + 0.1*(x1 + x2 + x3)
  
  # True outcome model
  y <- 2 + x1 + x2 + x3 + x4
  y <- add_noise(y,0.8)
  synthetic_pop <- as.data.frame(data.frame(y, x1, x2, x3, x4))
  colnames(synthetic_pop) <- c("y", "x1", "x2", "x3","x4")
  
  # Add ID column
  synthetic_pop$id <- 1:N
  
  return(synthetic_pop)
}

generate_synthetic_linear_high_noise <- function(N = 1000000) {
  z1 <- rbinom(N, 1, 0.5)
  z2 <- runif(N, 0, 2)
  z3 <- rexp(N, rate=1)
  z4 <- rchisq(N, df=4)

  x1 = z1
  x2 = z2 + 0.3*x1
  x3 = z3 + 0.2*(x1 + x2)
  x4 = z4 + 0.1*(x1 + x2 + x3)
  
  # True outcome model
  y <- 2 + x1 + x2 + x3 + x4 
  y <- add_noise(y,0.3)
  synthetic_pop <- as.data.frame(data.frame(y, x1, x2, x3, x4))
  colnames(synthetic_pop) <- c("y", "x1", "x2", "x3","x4")
  
  # Add ID column
  synthetic_pop$id <- 1:N
  
  return(synthetic_pop)
}

generate_synthetic_non_linear <- function(N = 1000000) {
  N = 1000000
  z1 <- rbinom(N, 1, 0.5)
  z2 <- runif(N, 0, 2)
  z3 <- rexp(N, rate=1)
  z4 <- rchisq(N, df=4)
  
  x1 = z1
  x2 = z2 + 0.3*x1
  x3 = z3 + 0.2*(x1 + x2)
  x4 = z4 + 0.1*(x1 + x2 + x3)
  
  # True outcome model
  raw_logit = (x2*x3 - (1/(x2*x4)))/(x1+1)
  
  raw_logit[raw_logit>quantile(raw_logit,0.99)] = quantile(raw_logit,0.99)
  raw_logit[raw_logit<quantile(raw_logit,0.01)] = quantile(raw_logit,0.01)
  
  y = sin(raw_logit)
  y <- add_noise(y,0.8)
  
  
  
  synthetic_pop <- as.data.frame(data.frame(y, x1, x2, x3, x4))
  colnames(synthetic_pop) <- c("y", "x1", "x2", "x3","x4")
  
  # Add ID column
  synthetic_pop$id <- 1:N
  
  return(synthetic_pop)
}

# POC Check model performance
# library(dplyr)
# library(xgboost)
# 
# kfolds = 5
# train_data = synthetic_pop %>% mutate(fold = sample(1:kfolds, size=nrow(synthetic_pop), replace=T))
# 
# cv_err_lm = list()
# cv_err_xgb = list()
# 
# i=1
# for(i in 1:kfolds){
#   in_data <- filter(train_data, fold!=i)
#   out_data <- filter(train_data, fold==i)
#   
#   pred_lm <- lm(y~x1+x2+x3+x4, data = in_data) %>% predict(out_data)
#   
#   dtrain <- xgb.DMatrix(data = as.matrix(in_data[,c("x1","x2","x3","x4")]), label= in_data$y)
#   model <- xgboost(data = dtrain, # the data   
#                    nround = 100, # max number of boosting iterations
#                    objective = "reg:squarederror",verbose = 0)  # the objective function
#   
#   pred_xgb = predict(model,as.matrix(out_data[,c("x1","x2","x3","x4")]))
#   
#   rmse_lm <- sqrt(mean((out_data$y - pred_lm)^2))
#   rmse_xgb <- sqrt(mean((out_data$y - pred_xgb)^2))
#   
#   r2_lm = cor(out_data$y,pred_lm)^2  
#   r2_xgb = cor(out_data$y,pred_xgb)^2  
#   
#   # Record the RMSE
#   cv_err_lm[[i]] <- list(r2 = r2_lm , rmse = rmse_lm)
#   cv_err_xgb[[i]] <- list(r2 = r2_xgb , rmse = rmse_xgb)
# }
# 
# do.call(rbind,cv_err_lm)
# do.call(rbind,cv_err_xgb)

# ===================================================================
# Generate and Save Datasets
# ===================================================================

# Create output directory if it doesn't exist
if(!dir.exists("data")) {
  dir.create("data")
}

# Generate synthetic dataset

list_dataset_name <- c("synthetic_linear_low_noise","synthetic_linear_high_noise","synthetic_non_linear")
cat("Generating synthetic dataset...\n")

for (dataset in list_dataset_name){
  generating_function <- get(paste0("generate_",dataset))
  synthetic_pop <- generating_function()
  
  save(synthetic_pop,file = paste0("data/",dataset,".RData"))
  write.csv(head(synthetic_pop, 100), file = paste0("data/",dataset,"_sample.csv"), row.names = FALSE)
  cat(paste0("\nSynthethic dataset saved to data/",dataset,".RData"))
}
cat("\nAll datasets generated successfully!\n")