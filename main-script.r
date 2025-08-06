# ===================================================================
# Main Script for Statistical Estimation Simulation
# ===================================================================
# This script contains all configurations and serves as the entry point
# for running the statistical estimation simulations.

# Load required libraries for parallelization and file management
library(parallel)
library(foreach)
library(doParallel)
library(pbapply)
library(dplyr)
library(tidyr)

# Load the core functions
source("code/statistical-estimation.R")
source("code/sampling-function.R")

# ===================================================================
# Configuration Setup
# ===================================================================

# Define dataset configurations
DATASETS <- list(
  synthethic_linear_low_noise = list(
    file_path = "data/synthetic_linear_low_noise.RData",
    data_object = "synthetic_pop",
    variables = list(
      outcome = "y",
      predictors = c("x1", "x2", "x3","x4")
    )
  ),
  synthethic_linear_high_noise = list(
    file_path = "data/synthetic_linear_high_noise.RData",
    data_object = "synthetic_pop",
    variables = list(
      outcome = "y",
      predictors = c("x1", "x2", "x3","x4")
    )
  ),
  synthethic_non_linear = list(
    file_path = "data/synthetic_non_linear.RData",
    data_object = "synthetic_pop",
    variables = list(
      outcome = "y",
      predictors = c("x1", "x2", "x3","x4")
    )
  )
)

# Define base configuration
BASE_CONFIG <- list(
  population_size = 1000000,     # Total sample size (N)
  ps_iterations = 20,            # Iterations for propensity score estimation
  random_seed = 42,              # Random seed for reproducibility
  replications = 150,             # Number of simulation replications
  output_dir = "result",        # Directory to save results
  run_parallel = T,
  n_cores = detectCores() - 1                    # Number of cores for parallel processing (fixed value)
)

# Create OKR configurations
all_config_okr <- expand.grid(
  dataset = "okr",
  sampling_technique = "okr",
  nonprob_size = c(1000,5000,25000),
  prob_size = c(500,1500,5000),
  stringsAsFactors = FALSE
)
all_config_okr$config_name <- paste0(all_config_okr$dataset, "_", 1:nrow(all_config_okr))
all_config_okr$description <- ""
list_configs_okr <- setNames(split(all_config_okr, seq(nrow(all_config_okr))), all_config_okr$config_name)

# Create synthetic data configurations
all_config_synthetic <- expand.grid(
  dataset = c("synthethic_linear_low_noise", "synthethic_linear_high_noise", "synthethic_non_linear"),
  sampling_technique = "synthetic",
  nonprob_size = c(1000,5000,25000),
  prob_size = c(500,1500,5000),
  method_variant = c("simple","complex"),
  stringsAsFactors = FALSE
)


all_config_synthetic$config_name <- paste0(
  all_config_synthetic$dataset, "_", "npsampling_", 
  all_config_synthetic$method_variant,"p",all_config_synthetic$prob_size,"np",all_config_synthetic$nonprob_size
)
all_config_synthetic$description <- ""
list_configs_syn <- setNames(split(all_config_synthetic, seq(nrow(all_config_synthetic))), all_config_synthetic$config_name)

# Define specific configurations to run
SIMULATION_CONFIGS <- c(list_configs_okr,list_configs_syn)
# Save configuration
save(SIMULATION_CONFIGS, file = "simulation_config_full.RData")

# ===================================================================
# File Management Functions
# ===================================================================

#' Create output directory if it doesn't exist
#' 
#' @param dir_path Path to directory
#' @return TRUE if directory exists or was created successfully
ensure_directory <- function(dir_path) {
  if(!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Created directory:", dir_path, "\n")
  }
  return(TRUE)
}

#' Generate a unique filename for results
#' 
#' @param config_name Name of the configuration
#' @param replication Replication number (optional)
#' @param output_dir Output directory
#' @return Full path to the output file
generate_filename <- function(config_name, replication = NULL, output_dir = "results") {
  if(is.null(replication)) {
    return(sprintf("%s/%s_combined.csv", output_dir, config_name))
  } else {
    return(sprintf("%s/%s_rep%03d.csv", output_dir, config_name, replication))
  }
}

#' Save simulation results to CSV file
#' 
#' @param results Results from a simulation
#' @param config_name Name of the configuration
#' @param config Configuration used
#' @param replication Replication number (optional)
#' @param output_dir Output directory
#' @return Path to the saved file
save_results <- function(results, config_name, config, replication = NULL, output_dir = "results") {
  # Ensure directory exists
  ensure_directory(output_dir)
  
  # Convert results to a data frame if it's a list
  if(is.list(results) && !is.data.frame(results)) {
    results_df <- do.call(rbind, results)
  } else {
    results_df <- results
  }
  
  # Add configuration info
  results_df$config <- config_name
  if(!is.null(replication)) {
    results_df$replication <- replication
  }
  
  # Add timestamp
  results_df$timestamp <- Sys.time()
  
  # Generate filename
  file_path <- generate_filename(config_name, replication, output_dir)
  
  # Save to CSV
  write.csv(results_df, file = file_path, row.names = FALSE)
  
  # Also save the configuration
  config_file <- sprintf("%s/%s_config.rds", output_dir, config_name)
  saveRDS(config, file = config_file)
  
  cat("Results saved to:", file_path, "\n")
  return(file_path)
}

#' Combine all replication results for a configuration
#' 
#' @param config_name Name of the configuration
#' @param max_replications Maximum number of replications
#' @param output_dir Output directory
#' @return Combined results data frame
combine_results <- function(config_name, max_replications, output_dir = "results") {
  all_results <- list()
  
  for(i in 1:max_replications) {
    file_path <- generate_filename(config_name, i, output_dir)
    if(file.exists(file_path)) {
      all_results[[i]] <- read.csv(file_path)
    }
  }
  
  if(length(all_results) > 0) {
    combined <- do.call(rbind, all_results)
    combined_file <- generate_filename(config_name, NULL, output_dir)
    write.csv(combined, file = combined_file, row.names = FALSE)
    cat("Combined results saved to:", combined_file, "\n")
    return(combined)
  } else {
    cat("No results found to combine for", config_name, "\n")
    return(NULL)
  }
}

# ===================================================================
# Simulation Runner Functions
# ===================================================================

#' Create a complete configuration by merging base config with specific config
#' 
#' @param base_config Base configuration
#' @param specific_config Specific configuration
#' @param config_name Name of the configuration
#' @return Complete configuration
create_complete_config <- function(base_config, specific_config, config_name) {
  # Start with base config
  config <- base_config
  
  # Add specific config settings
  for(name in names(specific_config)) {
    config[[name]] <- specific_config[[name]]
  }
  
  # Add dataset-specific info
  dataset_info <- DATASETS[[config$dataset]]
  if(is.null(dataset_info)) {
    stop(paste("Dataset not found:", config$dataset))
  }
  
  config$dataset_info <- dataset_info
  config$variables <- dataset_info$variables
  config$config_name <- config_name
  
  return(config)
}

#' Run a single replication of the simulation study
#'
#' @param config Configuration parameters
#' @return Data frame with estimation results
run_simulation_replication <- function(config) {
  # Extract configuration parameters
  outcome <- config$variables$outcome
  predictors <- config$variables$predictors
  
  # Load the specified dataset
  dataset_info <- DATASETS[[config$dataset]]
  load(dataset_info$file_path)
  
  # Get the population data object
  pop <- get(dataset_info$data_object)
  
  # Sample a subset of the population
  pop <- pop[sample(1:nrow(pop), size = config$population_size), ]
  if(!"id" %in% names(pop)) {
    pop$id <- 1:nrow(pop)
  }
  
  # Perform sampling using the dedicated sampling functions
  samples <- perform_sampling(pop, config)
  nonprob_sample <- samples$nonprob_sample
  prob_sample <- samples$prob_sample
  
  # Create outcome regression predictions
  or_predictions <- create_outcome_regression_predictions(
    nonprob_sample, prob_sample, predictors, outcome
  )
  
  # Calculate propensity scores with logistic regression
  nonprob_sample$ps_logreg_custom <- calculate_ps_logreg(
    nonprob_sample, prob_sample, predictors, pred_set = "np"
  )$ps_logreg_custom
  
  # Calculate propensity scores with XGBOOST
  nonprob_sample$ps_xgb_custom_1 <- calculate_ps_xgb_1(
    nonprob_sample, prob_sample, predictors, pred_set = "np"
  )$ps_xgb_custom
  
  nonprob_sample$ps_xgb_custom_2 <- calculate_ps_xgb_2(
    nonprob_sample, prob_sample, predictors, pred_set = "np"
  )$ps_xgb_custom
  
  # Cross validation for propensity score models
  cv_ps_logreg <- cross_validate_ps_model(nonprob_sample, prob_sample, predictors, calculate_ps_logreg)
  cv_ps_xgb_1 <- cross_validate_ps_model(nonprob_sample, prob_sample, predictors, calculate_ps_xgb_1)
  cv_ps_xgb_2 <- cross_validate_ps_model(nonprob_sample, prob_sample, predictors, calculate_ps_xgb_2)
  
  # Calculate pseudo-weight estimates (calibrated + uncalibrated)
  ps_calibrated_pw <- create_calibrated_pseudo_weight(
    nonprob_sample, prob_sample, predictors, pred_set = "np"
  )
  
  # Get pseudo-weight results
  pw_result <- pseudo_weight_estimates(nonprob_sample, ps_calibrated_pw, outcome)
  nonprob_sample = cbind(nonprob_sample,ps_calibrated_pw)
  
  # Cross validation for pseudo-weight methods
  cv_ps_cal_pw <- cross_validate_ps_model(nonprob_sample, prob_sample, predictors, create_calibrated_pseudo_weight)
  
  # Create all combinations of OR and PS models
  list_or_model <- names(or_predictions$nonprob_pred)
  list_ps_model <- c("ps_xgb_custom_1","ps_xgb_custom_2", "ps_logreg_custom")#,names(ps_calibrated_pw))
  combinations <- crossing(or = list_or_model, ps = list_ps_model)
  
  # Calculate double robust estimates
  all_result <- list()
  for (i in seq_len(nrow(combinations))) {
    all_result[[i]] <- create_dr_estimates(
      nonprob_sample = nonprob_sample,
      prob_sample = prob_sample,
      or_predictions = or_predictions,
      or_model = combinations$or[i],
      ps_model = combinations$ps[i],
      outcome = outcome,
      pop_size = config$population_size
    )
  }
  
  # Combine DR results
  dr_result <- combine_result_dr(all_result)
  
  # Combine all results
  # all_result <- rbind(dr_result, pw_result)
  all_result <- rbind(dr_result)
  all_result$y_true <- mean(samples$population[[outcome]])
  
  # Combine with metrics
  ps_metrics <- rbind(cv_ps_logreg, cv_ps_xgb_1,cv_ps_xgb_2,cv_ps_cal_pw)
  
  print(ps_metrics)
  all_result <- all_result %>%
  left_join(ps_metrics, by = join_by(ps_model == model)) %>%
  left_join(or_predictions$or_metrics, by = join_by(or_model == model))
  print(all_result)
  return(all_result)
}


# ===================================================================
# Run Specific Configuration
# ===================================================================

# Select which configurations to run (set to NULL to run all)
configs_to_run <- NULL  # Change this to run different configs

# If configs_to_run is NULL, run all configurations
if (is.null(configs_to_run)) {
  configs_to_run <- names(SIMULATION_CONFIGS)
}

# Filter the configurations to run
configs_subset <- SIMULATION_CONFIGS[names(SIMULATION_CONFIGS) %in% configs_to_run]

# Check if any valid configurations were found
if (length(configs_subset) == 0) {
  stop("No valid configurations selected. Check the config names.")
}

# Print summary of what will be run
cat("\n=======================================================\n")
cat("Statistical Estimation Simulation\n")
cat("=======================================================\n\n")
cat("Running the following configurations:\n")
for (config_name in names(configs_subset)) {
  config <- configs_subset[[config_name]]
  cat(sprintf("- %s: %s (dataset: %s, nonprob: %d, prob: %d)\n", 
              config_name, config$description, config$dataset, 
              config$nonprob_size, config$prob_size))
}
cat(sprintf("\nUsing %d cores for parallel processing\n", BASE_CONFIG$n_cores))
cat(sprintf("Results will be saved to '%s' directory\n", BASE_CONFIG$output_dir))

# Create log file
log_file <- paste0(BASE_CONFIG$output_dir, "/simulation_log.txt")
ensure_directory(BASE_CONFIG$output_dir)

# Set up cluster for parallel processing

if (BASE_CONFIG$run_parallel){
  n_cores <- BASE_CONFIG$n_cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
}

# Process each configuration
results <- list()

for(config_name in names(configs_subset)) {
  cat(paste("\n===== Processing configuration:", config_name, "=====\n"))
  
  # Create complete configuration
  config <- create_complete_config(BASE_CONFIG, configs_subset[[config_name]], config_name)
  # Setup output directory
  ensure_directory(config$output_dir)
  
  # Determine which replications still need to be run
  max_replications <- config$replications
  replications_to_run <- c()
  
  for(rep in 1:max_replications) {
    file_path <- generate_filename(config_name, rep, config$output_dir)
    if(!file.exists(file_path)) {
      replications_to_run <- c(replications_to_run, rep)
    }
  }
  
  if(length(replications_to_run) == 0) {
    cat(paste("All replications for", config_name, "already completed.\n"))
    
    # Combine results if needed
    combined_file <- generate_filename(config_name, NULL, config$output_dir)
    if(!file.exists(combined_file)) {
      combine_results(config_name, max_replications, config$output_dir)
    }
    
    results[[config_name]] <- read.csv(combined_file)
    next
  }
  
  cat(paste("Running", length(replications_to_run), "remaining replications for", config_name, "\n"))
  
  # Create a progress bar
  
  if (BASE_CONFIG$run_parallel){
    pb <- txtProgressBar(min = 0, max = length(replications_to_run), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Export required files and functions to the cluster
    clusterEvalQ(cl, {
      source("code/statistical-estimation.R")
      source("code/sampling-function.R")
      library(dplyr)
      library(tidyr)
      library(tidymodels)
      library(lightgbm)
      library(discrim)
      library(bonsai)
      library(yardstick)
    })
    
    # Export necessary objects to the cluster
    clusterExport(cl, c("DATASETS", "run_simulation_replication", "save_results"), envir = environment())
    
    # Run the replications in parallel
    foreach(rep = replications_to_run, .packages = c("dplyr", "tidyr"),.options.snow = opts) %dopar% {
      tryCatch({
        set.seed(config$random_seed + rep)
        result <- run_simulation_replication(config)
        save_results(result, config$config_name, config, rep, config$output_dir)
      }, error = function(e) {
        # Log the error
        error_message <- paste0("ERROR in replication ", rep, " for config ", config_name, ": ", conditionMessage(e))
        cat(error_message, "\n")
        # Write to a log file
        write(error_message, file = file.path(config$output_dir, "error_log.txt"), append = TRUE)
        # Return NULL to continue with next iteration
        return(NULL)
      })
    }
    
    close(pb)
  } else {
    for (rep in replications_to_run){
      set.seed(config$random_seed + rep)
      result <- run_simulation_replication(config)
      save_results(result, config$config_name, config, rep, config$output_dir)
      cat(paste("Finished with replication :",rep,config$config_name))
    }
  }

  # Combine results
  cat(paste("Combining results for", config_name, "\n"))
  combine_results(config_name, max_replications, config$output_dir)
  
  # Load the combined results
  combined_file <- generate_filename(config_name, NULL, config$output_dir)
  results[[config_name]] <- read.csv(combined_file)
}

cat("All configurations completed!\n")

cat("\n=======================================================\n")
cat("Simulation Completed\n")
cat("=======================================================\n\n")

# Stop the cluster
stopCluster(cl)


