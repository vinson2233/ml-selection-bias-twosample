# Correcting Selection Bias with Machine Learning

Author : [Vinson Ciawandy](https://www.linkedin.com/in/vinsonciawandy)  
Supervisors : 
- Prof. Ton de Waal (CBS)
- Dr. Sander Scholtus (CBS)
- Prof. Elise Dusseldorp (Leiden University)

This the code repository of the upcoming thesis master "Selection Bias Correction in a Two-Sample Setup Using Machine Learning".

The repository contains implementation statistical estimation using probabilistic and non-probability sampling methods, supporting various estimators including :
-  Double Robust (DR) estimators for nonprobability sample introduced by (Chen, 2022)[https://www.tandfonline.com/doi/abs/10.1080/01621459.2019.1677241]
-  Pseudo-weight estimators introduced by (Liu, 2023)[https://research.tilburguniversity.edu/en/publications/correcting-selection-bias-in-big-data-by-pseudo-weighting] 

Some novelty introduced with this work : 
- Applying Machine Learning algorithms with custom objective function to support Doubly Robust estimator 
- Applying Machine Learning algorithms on the Pseudo-Weight estimator

## Overview

The codebase developed with reproducibility, modularity, performance and maintanibility in mind.   
Features of this code base :
- Generate population's parameter estimation under two-sample setup (Non-probability sample and probability sampling)
- Easy integration of new dataset
- Customize data generation process in the simulation
- Flexibility to pick and choose which machine learning methods to use
- Multi-core support for monte carlo simulation
- Compare performance of different estimation methods

## Installation

### Requirements

- R 4.0.0 or later
- Required packages management handled through renv (see `Setup` section)

### Setup

1. Clone this repository:
   ```
   git clone https://gitea.cbsp.nl/VCNY/thesis_selection_bias_ml.git
   cd thesis_selection_bias_ml
   ```

2. Use renv to restore the exact package environment:
   ```R
   install.packages("renv")  # If you don't have renv installed yet
   renv::restore()
   ```
   This will automatically install all required packages with the exact versions used during development. There is a possibility that the library will ask whether you want to create a new project environment, or do you want to use the existing one. Just choose which one you prefer.

3. (Optional) Generate synthetic datasets:
   ```R
   source("sample-data-generator.R")
   ```
  
## File Structure
Inside the `code` folder, you can find : 
- `main-script.r`: Configuration and entry point to run simulations
- `statistical_estimation.R`: Core estimation functions and simulation logic
- `sampling_functions.R`: Dataset-specific sampling methods
- `sample-data-generator.R`: Script to generate synthetic test datasets if needed. Only need to be run once. It will produce `.Rdata` object.
- `results-analysis_4thJune2025.Rmd`: Analysis and visualization of simulation results.
- `renv.lock`: Package versions for reproducibility

![flow chart]("Diagram Study Flow.png")

## Usage

### Basic Usage

To run a simulation with default settings:

```R
source("main-script.R")
```

This will:
1. Load datasets specified in the configuration
2. Run specified number of replications
3. Generate estimates using various methods
4. Save results to the output directory

### Analyzing Results

To analyze simulation results:

```R
source("results-analysis.R")
```

This generates visualizations and summary statistics from the simulation results.

## Datasets and Sampling

The framework supports multiple datasets, each with their own sampling methods:

1. **OKR Dataset** 
   - Dutch Online Kilometer dataset from CBS

2. **Synthetic Dataset** (generated with `sample-data-generator.R`)
   - Synthetic data with known relationships.

### Adding New Datasets

To add a new dataset:

1. Add dataset configuration to `DATASETS` list in `main-script.R`:
   ```R
   DATASETS$new_dataset <- list(
     file_path = "data/new_dataset.RData",
     data_object = "new_dataset_pop",
     variables = list(
       outcome = "target_variable",
       predictors = c("pred1", "pred2", "pred3")
     )
   )
   ```

2. Create dataset-specific sampling functions in `sampling_functions.R`:
   ```R
   compute_prob_sample_probs_new_dataset <- function(data, target_size, predictors) {
     # Your probability sampling logic
   }
   
   compute_nonprob_sample_probs_new_dataset <- function(data, target_size, predictors) {
     # Your non-probability sampling logic
   }
   ```
The algorithm will find the suitable based on the function name, so make sure to make the dataset name consistent !

The framework is highly configurable through the `BASE_CONFIG` in `main-script.R`:

```R
BASE_CONFIG <- list(
  population_size = 1e5,           # Size of testing population
  sample_size = 50000,             # Total sample size
  ps_iterations = 20,              # Propensity score iterations
  replications = 100,              # Number of simulation replications
  output_dir = "results",          # Directory to save results
  n_cores = parallel::detectCores() - 1,  # Cores for parallel processing
  log_to_file = TRUE,              # Log output to file
  verbose_errors = TRUE,           # Detailed error information
)
```

### Simulation Configurations

Define specific scenarios to run:

```R
SIMULATION_CONFIGS <- list(
  config1 = list(
    dataset = "okr",
    nonprob_fraction = 0.2,
    prob_fraction = 0.2,
    description = "okr dataset with balanced sampling"
  ),
  config2 = list(
    dataset = "default",
    nonprob_fraction = 0.3,
    prob_fraction = 0.1,
    nonprob_method = "y",  # Use outcome-based sampling
    description = "Default dataset with outcome-based nonprob sampling"
  )
)
```

## Troubleshooting

### Error Logs
Detailed error logs are saved to:
```
results/logs/[config_name]/rep[XXX]_log.txt
```


<img src="https://www.rug.nl/digital-competence-centre/news/images/cbs-0.jpg" width="425"/> <img src="https://huisstijl.leidenuniv.nl/assets/files/ul-algemeen-internationaal-rgb-color.png" width="425"/> 


