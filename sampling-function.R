# ===================================================================
# Sampling Functions for Statistical Estimation
# ===================================================================
# This script contains functions for computing inclusion probabilities
# and performing sampling from a population, with dataset-specific logic.

library(dplyr)

# ===================================================================
# Generic Sampling Functions
# ===================================================================

#' Perform Poisson sampling from a population
#'
#' @param population Data frame containing the population
#' @param prob_col Column name containing inclusion probabilities
#' @return Data frame with the sampled observations
poisson_sampling <- function(population, prob_col) {
  inclusion_probs <- population[[prob_col]]
  included <- rbinom(length(inclusion_probs), 1, inclusion_probs)
  sampled_data <- population[included == 1, ]
  return(sampled_data)
}

#' Generic root-finding function for computing inclusion probabilities
#'s
#' @param data Data frame with population data
#' @param compute_pi Function to compute probabilities given constant c
#' @param target_size Target sample size
#' @return Constant c that gives the desired expected sample size
find_constant <- function(data, compute_pi, target_size) {
  cat("Find Constant Starting\n")
  
  # Define the intervals to try, from smallest to largest
  intervals <- list(
    c(-750,750),
    c(-350,350),
    c(-100,100)
  )
  
  result <- NULL
  success <- FALSE
  
  # Try each interval in sequence
  for (interval in intervals) {
    interval_text <- paste0("[", interval[1], ",", interval[2], "]")
    cat("Trying with interval", interval_text, "\n")
    
    tryCatch({
      result <- uniroot(
        function(const) sum(compute_pi(const, data)) - target_size,
        interval = interval,
        check.conv = TRUE,
        maxiter = 1000
      )$root
      
      # If we get here, the root finding was successful
      cat("Find Constant Success with interval", interval_text, "\n")

      success <- TRUE
      break  # Exit the loop once we find a solution
      
    }, error = function(e) {
      # Capture and report the error
      cat("Error with interval", interval_text, ":", conditionMessage(e), "\n")
      # No need to rethrow - we'll just try the next interval
    })
  }
  
  # Check if any interval was successful
  if (success) {
    return(result)
  } else {
    stop("Could not find root in any interval. Try specifying a wider range manually.")
  }
}

# Helper function to add noise
add_label_noise <- function(y, noise_ratio = 0.05, seed = NULL) {
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Ensure y is a vector
  y <- as.vector(y)
  
  # Create a copy of the original labels
  y_noisy <- y
  
  # Calculate number of labels to flip
  n_samples <- length(y)
  n_to_flip <- round(n_samples * noise_ratio)
  
  # Randomly select indices to flip
  indices_to_flip <- sample(n_samples, size = n_to_flip, replace = FALSE)
  
  # Flip the selected labels (assuming binary classification: 0 -> 1, 1 -> 0)
  y_noisy[indices_to_flip] <- 1 - y_noisy[indices_to_flip]
  
  # Return a list with the original features and the noisy labels
  return(y_noisy)
}
# ===================================================================
# Default Dataset (pop.RData) Sampling Functions
# ===================================================================

#' Compute inclusion probabilities for probability sample based on x1 (eerste_toel_dat_1)
#' for the default dataset
#'
#' @param data Data frame with population data
#' @param target_size Target sample size
#' @param predictor Name of the predictor variable to use
#' @return Data frame with added inclusion probabilities
compute_prob_sample_probs_okr <- function(data, target_size, predictor = "eerste_toel_dat_1") {
  # Calculate probabilities proportional to predictor
  compute_pi <- function(c2, df) {
    pi_p <- exp(c2 + df[[predictor]]/30) / (1 + exp(c2 + df[[predictor]]/30))
    return(pi_p)
  }
  
  # Find constant that gives desired expected sample size
  c2 <- find_constant(data, compute_pi, target_size)
  
  # Add inclusion probabilities to the data
  data$pi_p <- compute_pi(c2, data)
  data$meta_p_method <- "Prop_to_eerste_toel_dat_1"
  return(data)
}

#' Compute inclusion probabilities for non-probability sample based on x2-x3
#' for the default dataset
#'
#' @param data Data frame with population data
#' @param target_size Target sample size
#' @param predictors Vector of predictor variable names
#' @return Data frame with added inclusion probabilities
compute_nonprob_sample_probs_okr <- function(data, target_size, predictors = c("afl_gewkl", "leeftijd")) {
  # Calculate probabilities proportional to x2-x3/20
  compute_pi <- function(c1, df) {
    x <- df[[predictors[1]]] - df[[predictors[2]]] / 20
    pi_np <- exp(c1 + x) / (1 + exp(c1 + x))
    return(pi_np)
  }
  
  # Find constant that gives desired expected sample size
  c1 <- find_constant(data, compute_pi, target_size)
  
  # Add inclusion probabilities to the data
  data$pi_np <- compute_pi(c1, data)
  data$meta_np_method <- "Prop_to_afl_gewkl-leeftijd"
  return(data)
}


# ===================================================================
# Synthetic Dataset Sampling Functions
# ===================================================================

#' Compute inclusion probabilities for probability sample for synthetic dataset
#'
#' @param data Data frame with population data
#' @param target_size Target sample size
#' @param predictors Vector of predictor variable names
#' @return Data frame with added inclusion probabilities
compute_prob_sample_probs_synthetic <- function(data, target_size, predictors = c("x3")) {
  x3 = data$x3
  y = data$y
  
  # Calculate the inclusion probability pi_B (Probability sample)
  range_values = range(x3 + 0.03*y)
  min_val = range_values[1]
  max_val = range_values[2]
  c <- (max_val - 50*min_val)/49 # find c such that max(z)/min(z) = 50
  z = exp(c + x3 + 0.03*y)
  
  pi_B <- z / (1+z)
  pi_B = pi_B/sum(pi_B)*target_size # Normalize

  # Add inclusion probabilities to the data
  data$pi_p <- pi_B
  data$meta_np_method <- "proportional_to_x3+y"
  return(data)
}

#' Compute inclusion probabilities for non-probability sample for synthetic dataset
#'
#' @param data Data frame with population data
#' @param target_size Target sample size
#' @param predictors Vector of predictor variable names
#' @return Data frame with added inclusion probabilities
compute_nonprob_sample_probs_simple_synthetic <- function(data, target_size, predictors = c("x1","x2", "x3","x4")) {
  
  # Calculate the inclusion probability pi_A (Non-probability sample)
  compute_pi_A <- function(theta_0,data){
     
    x1 = data$x1 %>% scale()
    x2 = data$x2 %>% scale()
    x3 = data$x3 %>% scale()
    x4 = data$x4 %>% scale()
    
    logits = theta_0 + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4
    pi_A = 1 / (1 + exp(-logits))
    return(pi_A)
  }
  theta_0 = find_constant(data,compute_pi_A,target_size)

  # Add inclusion probabilities to the data
  data$pi_np <- compute_pi_A(theta_0,data)
  data$meta_np_method <- "simple"
  return(data)
}

#' Compute inclusion probabilities for non-probability sample for synthetic dataset
#'
#' @param data Data frame with population data
#' @param target_size Target sample size
#' @param predictors Vector of predictor variable names
#' @return Data frame with added inclusion probabilities
compute_nonprob_sample_probs_complex_synthetic <- function(data, target_size, predictors = c("x1","x2", "x3","x4")) {
  r = 50 # Ratio of maximum and minimum probability
  m = 1/r

  x1 = data$x1
  x2 = data$x2
  x3 = data$x3
  x4 = data$x4

  # True outcome model
  raw_logit = (x2*x3 - (1/(x2*x4)))/(x1+1)
  
  # Remove extreme outlier
  raw_logit[raw_logit>quantile(raw_logit,0.99)] = quantile(raw_logit,0.99)
  raw_logit[raw_logit<quantile(raw_logit,0.01)] = quantile(raw_logit,0.01)
  
  raw_logit = sqrt(abs(raw_logit)) # to spread the distribution
  logit = (raw_logit - min(raw_logit))/(max(raw_logit) - min(raw_logit))*pi/2 # Make the scale from 0 to pi/2
  
  y = sin(logit)
  scaled_pi_A = (y - min(y))/(max(y)-min(y))*(1-m) + m
  pi_A = target_size*(scaled_pi_A/sum(scaled_pi_A))
  
  # Add inclusion probabilities to the data
  data$pi_np <- pi_A
  data$meta_np_method <- "complex"
  return(data)
}


# ===================================================================
# Sampling Function Dispatcher
# ===================================================================

#' Get the appropriate sampling function for a dataset
#'
#' @param dataset_type Selection of sampling technique according to the dataset ("synthetic" or ""okr )
#' @param sample_type Type of sample ("prob" or "nonprob")
#' @param method_variant Optional variant of the method (e.g., "y" for outcome-based sampling)
#' @return Function to compute inclusion probabilities
get_sampling_function <- function(dataset_type, sample_type, method_variant = NULL) {
  # Function naming convention: compute_{type}_sample_probs_{dataset}_{variant}
  # where type is "prob" or "nonprob", dataset is the dataset name, and variant is optional
  
  # Build function name
  if(!is.null(method_variant)) {
    function_name <- paste0("compute_", sample_type, "_sample_probs_", method_variant, "_", dataset_type)
  } else {
    function_name <- paste0("compute_", sample_type, "_sample_probs_", dataset_type)
  }

  return(get(function_name))

}

#' Main sampling function that handles dataset-specific logic
#'
#' @param data Population data
#' @param config Configuration with dataset and sampling parameters
#' @return List with probability and non-probability samples
perform_sampling <- function(data, config) {
  # Extract parameters
  sampling_technique <- config$sampling_technique
  n_np <- config$nonprob_size
  n_p <- config$prob_size
  
  # Get sampling methods if specified
  prob_method <- if(!is.null(config$prob_method)) config$prob_method else NULL
  nonprob_method <- if(!is.null(config$nonprob_method)) config$nonprob_method else NULL
  method_variant <- if(!is.null(config$method_variant)) config$method_variant else NULL
  
  # Get appropriate sampling functions
  prob_sampling_func <- get_sampling_function(sampling_technique, "prob",)
  nonprob_sampling_func <- get_sampling_function(sampling_technique, "nonprob",method_variant)
  
  # Compute inclusion probabilities
  cat("Computing inclusion probabilities for probability sample...\n")
  data <- prob_sampling_func(data, n_p)
  
  cat("Computing inclusion probabilities for non-probability sample...\n")
  data <- nonprob_sampling_func(data, n_np)
  
  # Draw samples
  cat("Drawing probability sample...\n")
  prob_sample <- poisson_sampling(data, "pi_p")
  
  cat("Drawing non-probability sample...\n")
  nonprob_sample <- poisson_sampling(data, "pi_np")
  
  # Clean up: remove unnecessary columns from samples
  # nonprob_sample$pi_np <- NULL
  prob_sample$pi_np <- NULL
  
  return(list(
    prob_sample = prob_sample,
    nonprob_sample = nonprob_sample,
    population = data
  ))
}
