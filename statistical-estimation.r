# ===================================================================
# Statistical Estimation with Probabilistic Sampling
# ===================================================================
# Core implementation of statistical estimation methods

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidymodels)
library(probably)
library(xgboost)
library(discrim)
library(bonsai)
library(yardstick)
library(Rfast) # Fast matrix-vector operation

# Set consistent options
options(max.print = 1300)

# ===================================================================
# Outcome Regression (OR) Models
# ===================================================================

#' Create outcome regression predictions using various models
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Vector of predictor variable names
#' @param outcome Name of the outcome variable
#' @return List containing predictions for both samples
create_outcome_regression_predictions <- function(nonprob_sample, prob_sample, 
                                                 predictors, outcome) {
  
  # Create recipe
  or_recipe <- recipe(reformulate(predictors, outcome), data = nonprob_sample)
  
  # Define model specifications
  model_specs <- list(
    decision_tree = decision_tree() %>% set_engine("rpart") %>% set_mode("regression"),
    lin_reg = linear_reg() %>% set_engine("lm") %>% set_mode("regression"),
    random_forest = rand_forest(trees = 150) %>% set_engine("ranger") %>% set_mode("regression"),
    xgboost = boost_tree() %>% set_engine("xgboost") %>% set_mode("regression")
    
  )
  
  # Cross validation for metrics
  cv_folds <- vfold_cv(nonprob_sample, v = 5)
  
  # Function to fit a model with cross-validation and return metrics
  fit_with_cv <- function(model_spec, recipe) {
    workflow() %>% 
      add_recipe(recipe) %>% 
      add_model(model_spec) %>% 
      fit_resamples(
        resamples = cv_folds,
        metrics = metric_set(rmse, mae, rsq,msd),
        control = control_resamples(save_pred = TRUE)
      )
  }
  
  # Apply cross-validation to each model
  model_cv_results <- lapply(names(model_specs), function(model_name) {
    cv_result <- fit_with_cv(model_specs[[model_name]], or_recipe)
    metrics <- collect_metrics(cv_result)
    
    eval_result <- metrics %>% 
      mutate(model = model_name)
    
    return(eval_result)
  })
  
  # Combine results into a single dataframe
  cv_rmse_summary <- bind_rows(model_cv_results) %>% 
    dplyr::select(model, .metric, mean) %>%
    mutate(
      model = str_c("or_", model),
      .metric = str_c("or_", .metric)
    ) %>%
    pivot_wider(id_cols = model, names_from = .metric, values_from = mean)
  
  # Fit models
  model_fits <- lapply(model_specs, function(model_spec) {
    workflow() %>% 
      add_recipe(or_recipe) %>% 
      add_model(model_spec) %>% 
      fit(data = nonprob_sample)
  })
  
  # Make predictions for probability sample
  or_predictions_prob <- lapply(model_fits, predict, new_data = prob_sample)
  methods_name <- names(or_predictions_prob)
  or_predictions_prob <- as.data.frame(or_predictions_prob)
  colnames(or_predictions_prob) <- paste0("or_", methods_name)
  
  # Make predictions for non-probability sample
  or_predictions_nonprob <- lapply(model_fits, predict, new_data = nonprob_sample)
  methods_name <- names(or_predictions_nonprob)
  or_predictions_nonprob <- as.data.frame(or_predictions_nonprob)
  colnames(or_predictions_nonprob) <- paste0("or_", methods_name)
  
  return(list(
    nonprob_pred = or_predictions_nonprob,
    prob_pred = or_predictions_prob,
    or_metrics = cv_rmse_summary
  ))
}

# ===================================================================
# Propensity Score Models
# ===================================================================

#' Calculate metrics related to propensity model
#'
#' @param test_data Dataset with prediction inside of it
#' @param predictors Predictor variables
#' @return dataframe with metrics
propensity_metric_cal <- function(test_data, predictors) {
  NP <- filter(test_data, from_np == 1)
  P <- filter(test_data, from_np == 0)
  
  small_eps = 1e-5
  NP$pred[NP$pred==0] = small_eps
  P$pred[P$pred==1] = 1 - small_eps
  
  # Cross entropy
  if (nrow(P)!=0){
    entropy_np_term <- sum(log(NP$pred))
    entropy_p_term <- sum(log(1 - P$pred))
    n <- nrow(test_data)
    mxe <- -(1/n) * (entropy_np_term + entropy_p_term) 
    
    # Cal1
    w <- 1/NP$pred
    d <- 1/P$pi_p
    norm_constant <- sqrt((sapply(NP[predictors], var) + sapply(P[predictors], var))/2)
    cal1 <- sum((1/norm_constant) * abs(colSums(w*NP[, predictors])/sum(w) - 
                                          colSums(d*P[, predictors])/sum(d)))
  } else {
    mxe = -1
    cal1 = -1
  }
  
  # Brier score
  brier_score <- mean((test_data$pred - test_data$from_np)^2)
  
  return(data.frame(mxe = mxe, brier = brier_score, cal1 = cal1))
}

#' Run propensity score model cross validation
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Predictor variables
#' @param ps_function Propensity score function to evaluate
#' @return Cross validation result
cross_validate_ps_model <- function(nonprob_sample, prob_sample, predictors, ps_function) {
  nonprob_sample$intercept <- 1
  prob_sample$intercept <- 1
  k <- 5
  folds_np <- sample(1:k, nrow(nonprob_sample), replace = TRUE)
  folds_p <- sample(1:k, nrow(prob_sample), replace = TRUE)
  
  # Initialize performance metrics storage
  cv_results <- data.frame()
  for (i in 1:k) {
    # Split data
    train_nonprob_sample <- nonprob_sample[folds_np != i, ]
    train_prob_sample <- prob_sample[folds_p != i, ]
    
    test_data <- rbind(
      mutate(nonprob_sample[folds_np == i, c(predictors, "intercept","pi_p","id")], from_np = 1),
      mutate(prob_sample[folds_p == i, c(predictors, "intercept", "pi_p","id")], from_np = 0)
    )
    
    id_NP = test_data %>% filter(from_np==1) %>% .$id
    id_P = test_data %>% filter(from_np==0) %>% .$id
    id_intersect = intersect(id_NP,id_P)
    
    predictions <- ps_function(train_nonprob_sample, train_prob_sample, predictors, 
                              test_data, pred_set = "test")
    
    for (pred_col in names(predictions)) {
      
      test_data$pred <- predictions[[pred_col]]
      test_data_P_only = test_data %>% filter(id %in% id_P)
      test_data_no_intersection = test_data %>% filter(!id %in% id_intersect)
      
      # Evaluate metrics - 1
      ps_eval_result <- propensity_metric_cal(test_data_no_intersection, predictors)
      ps_eval_result$model <- pred_col
      ps_eval_result$subset <- "NP+P-Intersection"
      ps_eval_result$count_NP <- sum(test_data_no_intersection$from_np)
      ps_eval_result$count_P <- sum(1-test_data_no_intersection$from_np)
      ps_eval_result$count_intersection <- 0
      
      cv_results <- rbind(cv_results, ps_eval_result)
      
      # Evaluate metrics - 2
      ps_eval_result <- propensity_metric_cal(test_data_P_only, predictors)
      ps_eval_result$model <- pred_col
      ps_eval_result$subset <- "P"
      ps_eval_result$count_NP <- sum(test_data_P_only$from_np)
      ps_eval_result$count_P <- sum(1-test_data_P_only$from_np)
      ps_eval_result$count_intersection <- length(id_intersect)
      cv_results <- rbind(cv_results, ps_eval_result)
    }
  }
  
  cv_results <- cv_results %>% group_by(model,subset) %>% summarise_all(mean)
  return(cv_results) 
}

#' Fit propensity score model using logistic regression
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param aux_vars Variables to use for propensity modeling
#' @param max_iter Maximum number of iterations
#' @param plot_beta Whether to plot coefficient convergence
#' @return Fitted propensity score model
fit_propensity_score_lr <- function(nonprob_sample, prob_sample, aux_vars, 
                                   max_iter = 20, plot_beta = FALSE) {
  # Extract design matrices
  X_A <- as.matrix(nonprob_sample[, aux_vars])
  X_B <- as.matrix(prob_sample[, aux_vars])
  dB <- prob_sample$pi_p  # Inclusion probabilities
  max_iter=10
  # Initialize coefficients
  p <- ncol(X_A)
  beta <- as.matrix(rep(0, p), nrow = p)
  
  # Store beta history if plotting
  if (plot_beta) {
    betas <- beta
  }
  
  # Iterative estimation
  for (i in 1:max_iter) {
    # Calculate predicted probabilities
    pred <- 1 / (1 + exp(-(X_B %*% beta)))
    
    # Calculate gradient
    grad <- colSums(X_A) - colSums(matrix(rep(1 / dB * pred, p), ncol = p) * X_B)
    
    # Calculate Hessian
    W <- (1 / dB * pred * (1 - pred))
    hess <- t(X_B) %*% (matrix(rep(W, p), ncol = p) * X_B)
    
    # Update beta
    beta <- beta + solve(hess) %*% grad
    
    # Store beta if plotting
    if (plot_beta) {
      betas <- cbind(betas, beta)
    }
  }
  
  # Plot convergence if requested
  if (plot_beta) {
    plot <- t(betas) %>% 
      as.data.frame() %>% 
      mutate(iter = 1:(ncol(betas))) %>% 
      pivot_longer(cols = -iter) %>% 
      ggplot(aes(x = iter, y = value, col = name)) + 
      geom_line()
    print(plot)
  }
  
  # Create and return results
  result <- list(coef = beta, aux_vars = aux_vars, call = match.call())
  class(result) <- "ps_lr"
  return(result)
}

#' Predict method for propensity score logistic regression
#'
#' @param object Fitted propensity score model
#' @param newdata New data for prediction
#' @return Vector of propensity scores
predict.ps_lr <- function(object, newdata) {
  logit <- as.matrix(newdata[, object$aux_vars]) %*% object$coef
  ps <- 1 / (1 + exp(-logit))
  return(ps)
}

#' Calculate propensity scores using logistic regression
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Vector of predictor variable names
#' @param test_data Test data to be predicted if pred_set = "test"
#' @param pred_set One of c("test","np","np+p")
#' @param plot_beta Whether to plot coefficient convergence
#' @return Vector of propensity scores
calculate_ps_logreg <- function(nonprob_sample, prob_sample, predictors, 
                               test_data = NULL, pred_set = "test", plot_beta = FALSE) {
  nonprob_sample$intercept <- 1
  prob_sample$intercept <- 1
  
  model_ps_logreg <- fit_propensity_score_lr(
    nonprob_sample,
    prob_sample,
    aux_vars = c("intercept", predictors),
    plot_beta = plot_beta
  )
  
  if (pred_set == "np") {
    ps_logreg <- predict(model_ps_logreg, nonprob_sample)
  } else if (pred_set == "np+p") {
    ps_logreg <- predict(model_ps_logreg, rbind(
      nonprob_sample[c(predictors, "intercept")],
      prob_sample[c(predictors, "intercept")]
    ))
  } else if (pred_set == "test") {
    ps_logreg <- predict(model_ps_logreg, test_data)
  }
  
  return(data.frame("ps_logreg_custom" = ps_logreg))
}

#' Custom objective function for XGBoost propensity score model (version 1)
#'
#' @param preds Raw predictions
#' @param dtrain Training data
#' @return List with gradient and hessian
xgb_custom_obj_1 <- function(preds, dtrain) {
  # Transform raw predictions to probability
  p <- 1 / (1 + exp(-preds))

  # Create indicator for nonprob and prob samples
  IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
  d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample

  grad = - (IA - (1-IA)*d_B*p)
  hess = (1 - IA)*d_B*p*(1-p)

  return(list(grad = grad, hess = hess))
}

#' Custom objective function for XGBoost propensity score model (version 2)
#'
#' @param preds Raw predictions
#' @param dtrain Training data
#' @return List with gradient and hessian
xgb_custom_obj_2 <- function(preds, dtrain) {
  # Transform raw predictions to probability
  p <- 1 / (1 + exp(-preds))
  
  # Create indicator for nonprob and prob samples
  IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
  d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample
  
  # grad = - (IA - (1-IA)*d_B*p)
  # hess = (1 - IA)*d_B*p*(1-p)
  # # 
  # # hess[IA==1] <- 1e-6
  # 
  # # XGBoost expects gradient to be minimized
  
  grad <- ifelse(IA==1,
                 p - 1,      # decaying gradient for NP
                 d_B * p)      # original pseudo-likelihood for P
  hess <- ifelse(IA==1,
                 p * (1 - p),      # logistic hessian for NP
                 d_B * p * (1 - p))  # original for P
  return(list(grad = grad, hess = hess))
}



#' Calculate propensity scores using XGBoost (Version 1)
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Vector of predictor variable names
#' @param test_data Test data to be predicted if pred_set = "test"
#' @param pred_set One of c("test","np","np+p")
#' @return Vector of propensity scores
calculate_ps_xgb_1 <- function(nonprob_sample, prob_sample, predictors,
                             test_data = NULL, pred_set = "test") {
  X_all    <- as.matrix(rbind(nonprob_sample[predictors],
                              prob_sample[predictors]))
  y_all    <- c(rep(1, nrow(nonprob_sample)),
                rep(0, nrow(prob_sample)))
  w_all    <- c(rep(1, nrow(nonprob_sample)),
                1 / prob_sample$pi_p)
  
  # — set reproducible seed and compute split indices —
  n_total   <- nrow(X_all)
  train_pct <- 0.8
  n_train   <- floor(train_pct * n_total)
  
  # sample train rows
  train_idx <- sample(seq_len(n_total), size = n_train)
  val_idx   <- setdiff(seq_len(n_total), train_idx)
  
  # — build DMatrix objects —
  dtrain <- xgb.DMatrix(
    data   = X_all[train_idx, , drop = FALSE],
    label  = y_all[train_idx],
    weight = w_all[train_idx]
  )
  
  dval <- xgb.DMatrix(
    data   = X_all[val_idx, , drop = FALSE],
    label  = y_all[val_idx],
    weight = w_all[val_idx]
  )
  
  # — set up watchlist for early stopping etc. —
  watchlist <- list(
    train = dtrain,
    valid = dval
  )
  
  custom_eval <- function(preds, dtrain) {
    # Transform raw predictions to probability
    p <- 1 / (1 + exp(-preds))
    
    # Create indicator for nonprob and prob samples
    IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
    d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample
    
    value = sum(log(p[IA==1]/(1-p[IA==1]))) + sum(d_B[IA==0]*log(1-p[IA==0]))
    
    return(list(metric = "likelihood", value = value, higher_better = T))
  }
  
  
  # XGBoost parameters
  xgb_params <- list(
    max_depth = 3,
    min_child_weight = 5,
    subsample = 0.75,
    gamma = 1,
    eta = 0.03 # learning rate
  )
  
  # Train model with custom objective
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    early_stopping_rounds = 10,
    maximize = T,
    watchlist = watchlist,
    feval = custom_eval,
    nrounds = 250,
    obj = xgb_custom_obj_1,
    verbose = 0
  )
  
  # Make predictions based on requested prediction set
  if (pred_set == "np") {
    pred_data <- as.matrix(nonprob_sample[predictors])
  } else if (pred_set == "np+p") {
    pred_data <- as.matrix(rbind(
      nonprob_sample[predictors],
      prob_sample[predictors]
    ))
  } else if (pred_set == "test") {
    pred_data <- as.matrix(test_data[predictors])
  }
  
  # To get raw scores, use output_margin = TRUE
  pred <- predict(model, newdata = xgb.DMatrix(pred_data), outputmargin = F)
  
  # Convert to probabilities
  ps_xgb <- 1 / (1 + exp(-pred))
  return(data.frame(ps_xgb_custom_1 = ps_xgb))
}

#' Calculate propensity scores using XGBoost (Version 2)
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Vector of predictor variable names
#' @param test_data Test data to be predicted if pred_set = "test"
#' @param pred_set One of c("test","np","np+p")
#' @return Vector of propensity scores
calculate_ps_xgb_2 <- function(nonprob_sample, prob_sample, predictors,
                               test_data = NULL, pred_set = "test") {
  X_all    <- as.matrix(rbind(nonprob_sample[predictors],
                              prob_sample[predictors]))
  y_all    <- c(rep(1, nrow(nonprob_sample)),
                rep(0, nrow(prob_sample)))
  w_all    <- c(rep(1, nrow(nonprob_sample)),
                1 / prob_sample$pi_p)
  
  # — set reproducible seed and compute split indices —
  n_total   <- nrow(X_all)
  train_pct <- 0.8
  n_train   <- floor(train_pct * n_total)
  
  # sample train rows
  train_idx <- sample(seq_len(n_total), size = n_train)
  val_idx   <- setdiff(seq_len(n_total), train_idx)
  
  # — build DMatrix objects —
  dtrain <- xgb.DMatrix(
    data   = X_all[train_idx, , drop = FALSE],
    label  = y_all[train_idx],
    weight = w_all[train_idx]
  )
  
  dval <- xgb.DMatrix(
    data   = X_all[val_idx, , drop = FALSE],
    label  = y_all[val_idx],
    weight = w_all[val_idx]
  )
  
  # — set up watchlist for early stopping etc. —
  watchlist <- list(
    train = dtrain,
    valid = dval
  )
  
  custom_eval <- function(preds, dtrain) {
    # Transform raw predictions to probability
    p <- 1 / (1 + exp(-preds))
    
    # Create indicator for nonprob and prob samples
    IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
    d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample
    
    value = sum(log(p[IA==1]/(1-p[IA==1]))) + sum(d_B[IA==0]*log(1-p[IA==0]))
    
    return(list(metric = "likelihood", value = value, higher_better = T))
  }
  
  
  # XGBoost parameters
  xgb_params <- list(
    max_depth = 3,
    min_child_weight = 5,
    subsample = 0.75,
    gamma = 1,
    eta = 0.03 # learning rate
  )
  
  # Train model with custom objective
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    early_stopping_rounds = 10,
    maximize = T,
    watchlist = watchlist,
    feval = custom_eval,
    nrounds = 250,
    obj = xgb_custom_obj_2,
    verbose = 0
  )
  
  # Make predictions based on requested prediction set
  if (pred_set == "np") {
    pred_data <- as.matrix(nonprob_sample[predictors])
  } else if (pred_set == "np+p") {
    pred_data <- as.matrix(rbind(
      nonprob_sample[predictors],
      prob_sample[predictors]
    ))
  } else if (pred_set == "test") {
    pred_data <- as.matrix(test_data[predictors])
  }
  
  # To get raw scores, use output_margin = TRUE
  pred <- predict(model, newdata = xgb.DMatrix(pred_data), outputmargin = F)
  
  # Convert to probabilities
  ps_xgb <- 1 / (1 + exp(-pred))
  return(data.frame(ps_xgb_custom_2 = ps_xgb))
}



#' Custom objective function for XGBoost propensity score model (version 3)
#'
#' @param preds Raw predictions
#' @param dtrain Training data
#' @return List with gradient and hessian
xgb_custom_obj_3 <- function(preds, dtrain) {
  # Transform raw predictions to probability
  p <- 1 / (1 + exp(-preds))
  
  # Create indicator for nonprob and prob samples
  IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
  d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample
  
  
  grad <- ifelse(IA==1,
                 d_B * (p - 1),      # decaying gradient for NP
                 d_B * (p - 0))      # original pseudo-likelihood for P
  hess <- ifelse(IA==1,
                 d_B * p * (1 - p),  # logistic hessian for NP
                 d_B * p * (1 - p))  # original for P
  return(list(grad = grad, hess = hess))
}

#' Calculate propensity scores using XGBoost (Version 3)
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param predictors Vector of predictor variable names
#' @param test_data Test data to be predicted if pred_set = "test"
#' @param pred_set One of c("test","np","np+p")
#' @return Vector of propensity scores
calculate_ps_xgb_3 <- function(nonprob_sample, prob_sample, predictors,
                               test_data = NULL, pred_set = "test") {
  X_all    <- as.matrix(rbind(nonprob_sample[predictors],
                              prob_sample[predictors]))
  y_all    <- c(rep(1, nrow(nonprob_sample)),
                rep(0, nrow(prob_sample)))
  # w_all    <- c(rep(1, nrow(nonprob_sample)),
                # 1 / prob_sample$pi_p)
  ps_xgb = rep(1,nrow(nonprob_sample))
  
  
  # — set reproducible seed and compute split indices —
  n_total   <- nrow(X_all)
  train_pct <- 0.8
  n_train   <- floor(train_pct * n_total)
  
  # sample train rows
  train_idx <- sample(seq_len(n_total), size = n_train)
  val_idx   <- setdiff(seq_len(n_total), train_idx)
  
  df_result = data.frame()
  for( i in 1:100){
  print(i)
  w_all    <- c(1/ps_xgb,
                1 / prob_sample$pi_p)
  
  # — build DMatrix objects —
  dtrain <- xgb.DMatrix(
    data   = X_all[train_idx, , drop = FALSE],
    label  = y_all[train_idx],
    weight = w_all[train_idx]
  )
  
  dval <- xgb.DMatrix(
    data   = X_all[val_idx, , drop = FALSE],
    label  = y_all[val_idx],
    weight = w_all[val_idx]
  )
  
  # — set up watchlist for early stopping etc. —
  watchlist <- list(
    train = dtrain,
    valid = dval
  )
  
  custom_eval <- function(preds, dtrain) {
    # Transform raw predictions to probability
    p <- 1 / (1 + exp(-preds))
    
    # Create indicator for nonprob and prob samples
    IA <-  getinfo(dtrain, "label")  # 1 for nonprob sample, 0 for prob sample
    d_B <- getinfo(dtrain, "weight")  # Inclusion weights for prob sample
    
    value = sum(log(p[IA==1]/(1-p[IA==1]))) + sum(d_B[IA==0]*log(1-p[IA==0]))
    
    return(list(metric = "likelihood", value = value, higher_better = T))
  }
  
  
  # XGBoost parameters
  xgb_params <- list(
    max_depth = 3,
    min_child_weight = 5,
    subsample = 0.75,
    gamma = 1,
    eta = 0.03 # learning rate
  )
  
  # Train model with custom objective
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    early_stopping_rounds = 10,
    maximize = T,
    watchlist = watchlist,
    feval = custom_eval,
    nrounds = 250,
    obj = xgb_custom_obj_3,
    verbose = 1
  )
  
  # Make predictions based on requested prediction set
  if (pred_set == "np") {
    pred_data <- as.matrix(nonprob_sample[predictors])
  } else if (pred_set == "np+p") {
    pred_data <- as.matrix(rbind(
      nonprob_sample[predictors],
      prob_sample[predictors]
    ))
  } else if (pred_set == "test") {
    pred_data <- as.matrix(test_data[predictors])
  }
  
  # To get raw scores, use output_margin = TRUE
  pred <- predict(model, newdata = xgb.DMatrix(pred_data), outputmargin = F)
  
  # Convert to probabilities
  ps_xgb <- 1 / (1 + exp(-pred))
  
  temp = data.frame(iteration = i,id = nonprob_sample$id, pi_np_hat = ps_xgb)
  df_result = rbind(df_result,temp)
  }
  return(data.frame(ps_xgb_custom = ps_xgb))
}

# idx = nonprob_sample$id %>% sample(16)
# true_prob = nonprob_sample %>% filter(id %in% idx ) %>% select(id,pi_np)
# df_result %>% left_join(true_prob) %>% filter(train_idx%in%idx) %>% ggplot(aes(x=iteration,y=pi_np_hat)) + geom_line() + facet_wrap(~id)+ geom_hline(aes(yintercept = pi_np))
# 
# nonprob_sample$pi_np_hat = ps_xgb
# ggplot(nonprob_sample,aes(x=pi_np,y=pi_np_hat)) + geom_point(alpha=0.1) + geom_abline(slope=1,intercept = 0)

# ===================================================================
# Pseudo-Weight Methods
# ===================================================================

create_calibrated_pseudo_weight <- function(
    nonprob_sample, prob_sample, predictors,
    test_data = NULL, pred_set = c("test", "np"),
    v = 5,
    do_platt   = TRUE,
    do_isotonic = TRUE,
    smooth_platt = FALSE        # TRUE = GAM–smoothed Platt
) {
  
  stopifnot(do_platt | do_isotonic)   # must pick at least one method
  pred_set <- match.arg(pred_set)
  
  ## ---------------------------------------------------------------- ##
  ## 0.  Labelled data                                                ##
  ## ---------------------------------------------------------------- ##
  nonprob_sample$S_star <- 1
  prob_sample$S_star    <- 0
  
  combined <- dplyr::bind_rows(
    nonprob_sample[c(predictors, "id", "S_star")],
    prob_sample[c(predictors, "id", "S_star")]
  ) |>
    dplyr::anti_join(
      tibble(id = intersect(nonprob_sample$id, prob_sample$id)),
      by = "id"
    ) |>
    mutate(S_star = factor(S_star))
  
  ## ---------------------------------------------------------------- ##
  ## 1.  v‑fold CV                                                    ##
  ## ---------------------------------------------------------------- ##
  folds <- vfold_cv(combined, v = v, strata = S_star)
  
  ## ---------------------------------------------------------------- ##
  ## 2.  Common recipe                                                ##
  ## ---------------------------------------------------------------- ##
  base_recipe <- recipe(reformulate(predictors,response="S_star"), data = combined) #|>
  #    update_role(id, new_role = "id variable")
  
  ## ---------------------------------------------------------------- ##
  ## 3.  Model specs                                                  ##
  ## ---------------------------------------------------------------- ##
  model_specs <- list(
    logreg         = logistic_reg()              |> set_engine("glm"),
    qda_model      = discrim_quad()            |> set_engine("MASS"),
    decision_tree  = decision_tree()             |> set_engine("rpart"),
    random_forest  = rand_forest(trees = 150)    |> set_engine("ranger"),
    xgboost           = boost_tree()                |> set_engine("xgboost")
  ) |>
    purrr::map(~ .x |> set_mode("classification"))
  
  ## helper to build workflow ---------------------------------------- ##
  wflow <- \(spec) workflow() |> add_recipe(base_recipe) |> add_model(spec)
  
  ## target data to predict on --------------------------------------- ##
  pred_data <- if (pred_set == "np") nonprob_sample[c(predictors,"pi_p")] else test_data
  
  
  ## ---------------------------------------------------------------- ##
  ## 4.  For each algorithm …                                         ##
  ## ---------------------------------------------------------------- ##
  out <- purrr::imap_dfc(
    model_specs,
    function(spec, name) {
      ## 4‑a  OOF predictions via fit_resamples ---------------------- ##
      cv_fit <- fit_resamples(
        wflow(spec), resamples = folds,
        control = control_resamples(save_pred = TRUE)
      )
      oof <- collect_predictions(cv_fit, summarize = FALSE)
      
      ## 4‑b  Calibrators ------------------------------------------- ##
      cal_objs <- list()
      
      if (do_platt) {
        cal_objs$platt <- cal_estimate_logistic(oof,S_star,
                                                smooth   = smooth_platt
        )
      }
      if (do_isotonic) {
        cal_objs$iso <- cal_estimate_isotonic(oof,S_star
        )
      }
      
      ## 4‑c  Fit final model on the full data ----------------------- ##
      final_fit <- fit(wflow(spec), combined)
      
      raw <- predict(final_fit, pred_data, type = "prob") #|>
      #        transmute(!!paste0("uncal_", name) := .pred_1)
      
      ## apply each calibrator -------------------------------------- ##
      cal_cols <- purrr::imap_dfc(
        cal_objs,
        \(obj, lab) {
          cal_apply(raw,obj) |> transmute(!!paste0("cal_", lab, "_", name) := .pred_1) 
        }
      )
      
      bind_cols(raw %>% transmute(!!paste0("uncal_",name) := .pred_1), cal_cols)
    }
  )
  
  ps_predictions = as.data.frame(out)
  
  limit = 0.1
  ps_predictions[ps_predictions==1] = 1 - limit
  ps_predictions[ps_predictions==0] = limit
  
  O_hat <- ps_predictions / (1 - ps_predictions)
  pi_np_hat = O_hat / (O_hat + 1/pred_data$pi_p - 1)
  return(pi_np_hat)
}

# ===================================================================
# Estimation Methods
# ===================================================================

#' Create double robust estimates using various combinations of models
#'
#' @param nonprob_sample Non-probability sample data
#' @param prob_sample Probability sample data
#' @param or_predictions Outcome regression predictions
#' @param or_model Name of the OR model to use
#' @param ps_model Name of the PS model column in nonprob_sample
#' @param outcome Name of the outcome variable
#' @param pop_size Total population size
#' @return Data frame with estimation results
create_dr_estimates <- function(nonprob_sample, prob_sample, or_predictions,
                               or_model, ps_model, outcome, pop_size) {
  # Extract predictions and weights
  m_SA <- or_predictions$nonprob_pred[[or_model]]
  m_SB <- or_predictions$prob_pred[[or_model]]
  d_A <- 1 / nonprob_sample[[ps_model]]
  d_B <- 1 / prob_sample$pi_p
  
  # Calculate population size estimates
  N_A <- sum(d_A)
  N_B <- sum(d_B)
  N <- pop_size
  
  # Calculate various estimators
  mu_A <- mean(nonprob_sample[[outcome]])  # Naive estimator
  mu_reg <- (1 / N_B) * sum(d_B * m_SB)  # Regression estimator
  
  # Inverse probability weighting estimators
  mu_IPW1 <- (1 / N) * sum(nonprob_sample[[outcome]] * d_A)
  mu_IPW2 <- (1 / N_A) * sum(nonprob_sample[[outcome]] * d_A)
  
  # Double robust estimators
  mu_DR1 <- (1 / N) * sum(d_A * (nonprob_sample[[outcome]] - m_SA)) + (1 / N) * sum(d_B * m_SB)
  mu_DR2 <- (1 / N_A) * sum(d_A * (nonprob_sample[[outcome]] - m_SA)) + (1 / N_B) * sum(d_B * m_SB)
  
  # Combine results
  df_result <- rbind(
    c("-", "-", "naive", mu_A),
    c("-", or_model, "REG", mu_reg),
    c(ps_model, "-", "IPW1", mu_IPW1),
    c(ps_model, "-", "IPW2", mu_IPW2),
    c(ps_model, or_model, "DR1", mu_DR1),
    c(ps_model, or_model, "DR2", mu_DR2)
  ) %>% data.frame()
  
  # Add column names
  colnames(df_result) <- c("ps_model", "or_model", "method", "estimate")
  return(df_result)
}

#' Combine results from multiple double robust estimations
#'
#' @param all_results List of result data frames
#' @return Combined data frame of results
combine_result_dr <- function(all_results) {
  all_results <- do.call(rbind, all_results)
  all_results$estimate <- as.numeric(all_results$estimate)
  all_results <- all_results[!duplicated(all_results), ]
  return(all_results)
}

#' Calculate pseudo-weight estimates for multiple propensity models
#'
#' @param nonprob_sample Non-probability sample data
#' @param ps_predictions Data frame with propensity score predictions
#' @param outcome Name of the outcome variable
#' @return Data frame with estimation results
pseudo_weight_estimates <- function(nonprob_sample, ps_predictions, outcome) {
  # # Calculate odds ratios
  list_methods <- names(ps_predictions)
  # 
  # Initialize results data frame
  result <- data.frame()
  
  # Calculate estimates for each method
  for (method in list_methods) {
    # Calculate weights
    w <- 1/ps_predictions[[method]]
    # Calculate weighted mean
    mu_pw <- sum(w * nonprob_sample[[outcome]]) / sum(w)
    
    # Add to results
    result <- rbind(result, c(method, "-", "pseudo-weighting", mu_pw))
  }
  
  # Add column names and convert estimate to numeric
  names(result) <- c("ps_model", "or_model", "method", "estimate")
  result$estimate <- as.numeric(result$estimate)
  
  return(result)
}
