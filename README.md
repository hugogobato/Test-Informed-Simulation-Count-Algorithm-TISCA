# TISCA: Test-Informed Simulation Count Algorithm

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

The Test-Informed Simulation Count Algorithm (TISCA) is a statistically rigorous framework for determining the optimal number of simulation replications needed in research studies. Instead of selecting an arbitrary number of simulations, TISCA iteratively runs simulations and performs statistical tests until a pre-specified power is achieved for detecting user-defined minimum detectable effect sizes.

TISCA addresses a critical gap in simulation-based research: ensuring that studies have sufficient statistical power while avoiding computational waste. By dynamically determining the required number of simulations, TISCA enhances the reliability and reproducibility of research findings.

## Key Features

- **Statistical Rigor**: Ensures simulations have sufficient power to detect meaningful differences
- **Resource Efficiency**: Prevents wasteful computation by running only the necessary number of simulations
- **Reproducibility**: Standardizes the approach to determining simulation counts
- **Comprehensive Reporting**: Generates detailed statistics and visualizations of power and p-values
- **Multiple Testing Support**: Includes various methods for p-value adjustment
- **Flexible Implementation**: Works with any simulation function that returns performance metrics

## Installation

### R Implementation

```R
# Clone this repository
git clone https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git

# Source the TISCA.R file in your R script
source("path/to/TISCA.R")
```

### Python Implementation (Beta)

A Python implementation of TISCA is also available, though it is currently in beta phase:

```python
# Clone this repository
git clone https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git

# Import the TISCA module in your Python script
from tisca import run_tisca
```

**Note:** The Python implementation is still in beta testing and may have limitations compared to the R version. Additionally, the fine-tuned LLM available on the Streamlit website currently only supports the R implementation.

## Required Packages

### For R Implementation

TISCA depends on the following R packages:

```R
# Install required packages
install.packages(c("stats", "utils", "progress", "ggplot2"))
```

### For Python Implementation

The Python version requires the following packages:

```python
# Install required packages
pip install numpy pandas scipy tqdm matplotlib seaborn statsmodels
```

## Basic Usage

### R Implementation

Here's a simple example of how to use TISCA in R:

```R
# Source the TISCA.R file
source("TISCA.R")

# Define a simulation function that returns performance metrics
sim_func <- function(seed) {
  set.seed(seed)
  # Run your simulation and return metrics as a data frame
  data.frame(
    proposed_metric = rnorm(1, 0.5, 0.1),
    benchmark_metric = rnorm(1, 0.6, 0.1)
  )
}

# Define comparison pairs (which metrics to compare)
comparison_pairs <- list(
  c("proposed_metric", "benchmark_metric")
)

# Define minimum detectable effect sizes
# Use negative values for metrics where lower is better (e.g., RMSE)
# Use positive values for metrics where higher is better (e.g., coverage)
mdes <- c(-0.1)  # Detect a difference of at least 0.1 in favor of proposed_metric

# Run TISCA
results <- run_tisca(
  sim_func = sim_func,
  comparison_pairs = comparison_pairs,
  mdes = mdes,
  target_power = 0.80,
  alpha = 0.05,
  batch_size = 10,
  initial_count = 10,
  correction_method = "none",
  verbose = TRUE,
  save_results = TRUE
)

# Access results
print(results$J_final)  # Final number of simulations required
print(results$summary_table)  # Summary of results
```

### Python Implementation

See in this Github Repository in folder Usage_Examples/Synthetic_Data_Example_Python/ .

## Parameters

The `run_tisca` function accepts the following parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sim_func` | Function that runs one simulation and returns metrics | (Required) |
| `comparison_pairs` | List of metric pairs to compare | (Required) |
| `mdes` | Minimum detectable effect sizes | (Required) |
| `target_power` | Desired statistical power | 0.80 |
| `alpha` | Significance level | 0.05 |
| `batch_size` | Number of simulations per batch | 50 |
| `initial_count` | Initial simulations before first power check | batch_size |
| `correction_method` | Method for p-value adjustment | "none" |
| `verbose` | Whether to print progress updates | TRUE |
| `save_results` | Whether to save results to CSV files | TRUE |
| `output_dir` | Directory to save output files | "." |

## Return Value

The `run_tisca` function returns a list containing:

- `J_final`: The final number of simulations required
- `P_raw`: Vector of raw p-values
- `P_adj`: Vector of adjusted p-values
- `Stats`: Vector of test statistics
- `P_achieved`: Vector of achieved powers
- `results`: Data frame with all simulation results
- `comparison_pairs`: The comparison pairs used
- `mdes`: The minimum detectable effect sizes used
- `summary_table`: A data frame summarizing the results
- `power_tracking`: Data frame tracking power estimates across simulation counts
- `pvalue_tracking`: Data frame tracking p-values across simulation counts

## Advanced Example: Treatment Effect Estimation

TISCA is particularly useful for comparing methods in causal inference and treatment effect estimation. Here's an example comparing different methods for estimating heterogeneous treatment effects:

```R
# Define a simulation function for treatment effect estimation
simulate_treatment_effects <- function(seed) {
  set.seed(seed)
  
  # Generate synthetic data
  n <- 100
  X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  tau <- 2 * X[,1] + X[,2]^2  # True treatment effects
  
  # Fit models (proposed and benchmark)
  proposed_results <- fit_proposed_model(X, tau)
  benchmark_results <- fit_benchmark_model(X, tau)
  
  # Calculate performance metrics
  proposed_pehe <- sqrt(mean((proposed_results$point_estimates - tau)^2))
  benchmark_pehe <- sqrt(mean((benchmark_results$point_estimates - tau)^2))
  
  proposed_coverage <- calculate_coverage(proposed_results$posterior_samples, tau)
  benchmark_coverage <- calculate_coverage(benchmark_results$posterior_samples, tau)
  
  # Return results as a data frame
  return(data.frame(
    proposed_pehe = proposed_pehe,
    benchmark_pehe = benchmark_pehe,
    proposed_coverage = proposed_coverage,
    benchmark_coverage = benchmark_coverage
  ))
}

# Define comparison pairs
comparison_pairs <- list(
  c("proposed_pehe", "benchmark_pehe"),         # Compare PEHE (lower is better)
  c("proposed_coverage", "benchmark_coverage")  # Compare coverage (higher is better)
)

# Define minimum detectable effect sizes
mdes <- c(-0.1, 0.05)  # -0.1 for PEHE (lower is better), 0.05 for coverage (higher is better)

# Run TISCA
results <- run_tisca(
  sim_func = simulate_treatment_effects,
  comparison_pairs = comparison_pairs,
  mdes = mdes,
  target_power = 0.90,
  alpha = 0.01,
  batch_size = 10,
  correction_method = "holm"
)
```

## Visualizations

TISCA generates two main visualizations when `save_results = TRUE`:

1. **Power vs. Number of Simulations**: Shows how the estimated power for each comparison changes as more simulations are run.

2. **Significance of Comparisons**: A bar plot showing the significance (-log10 p-value) of each comparison after the final number of simulations.

These visualizations help researchers understand the convergence of power estimates and the statistical significance of their findings.

## Implementation Tips

1. **Simulation Function Design**: Ensure your simulation function accepts a seed parameter and returns a data frame with all metrics for comparison.

2. **Choosing MDEs**: Select minimum detectable effect sizes based on what would be practically significant in your field.

3. **Batch Size**: For computationally expensive simulations, use smaller batch sizes to avoid long waiting times between power checks.

4. **Multiple Testing**: When comparing many metrics, consider using a correction method like "holm" or "BH" to control the family-wise error rate or false discovery rate.

5. **Parallel Processing**: For very expensive simulations, consider implementing parallel processing within your simulation function.

## TISCA Streamlit App

We've created a companion web application to help researchers implement TISCA in their simulation studies. Visit [https://tisca-llm-app.streamlit.app/](https://tisca-llm-app.streamlit.app/) to access a fine-tuned LLM that can assist with:

- Structuring your R simulation code to be compatible with TISCA's `sim_func` input
- Guidance on choosing appropriate TISCA parameters (MDEs, comparison_pairs, alpha, target_power)
- Helping you interpret TISCA outputs
- Assisting with debugging R code related to your TISCA implementation

The web application provides an interactive interface for researchers who may be new to TISCA or need assistance with implementation details.

**Note:** Currently, the fine-tuned LLM only supports the R implementation of TISCA. Support for the Python implementation is planned for future updates.

## Citation

If you use TISCA in your research, please cite:

```
@misc{https://doi.org/10.48550/arxiv.2409.05161,
  doi = {10.48550/ARXIV.2409.05161},
  url = {https://arxiv.org/abs/2409.05161},
  author = {Souto,  Hugo Gobato and Neto,  Francisco Louzada},
  keywords = {Methodology (stat.ME),  FOS: Computer and information sciences,  FOS: Computer and information sciences},
  title = {Beyond Arbitrary Replications: A Principled Approach to Simulation Design in Causal Inference},
  publisher = {arXiv},
  year = {2024},
  copyright = {Creative Commons Attribution Non Commercial Share Alike 4.0 International}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- This work was inspired by the need for more rigorous simulation count determination in causal inference research
- Special thanks to all contributors and users who have provided valuable feedback


