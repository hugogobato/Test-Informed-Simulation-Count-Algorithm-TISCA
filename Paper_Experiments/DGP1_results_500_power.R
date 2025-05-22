required_packages <- c("car",'RVAideMemoire')

# Function to install and load packages
install_load_packages <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
    library(package, character.only = TRUE)
  }
}

# Install and load the packages
install_load_packages(required_packages)


#' Calculate Power for Welch's t-test via Simulation
#'
#' Estimates the power of Welch's two-sample t-test (unequal variances)
#' for detecting a specific true difference in means (`delta`), based on the
#' characteristics (sample sizes and standard deviations) of the provided data groups.
#'
#' @param group1_values Numeric vector of observed data for group 1.
#' @param group2_values Numeric vector of observed data for group 2.
#' @param delta The true difference in means (group2_mean - group1_mean) under the
#'   alternative hypothesis that you want to detect.
#' @param sig.level The significance level (alpha) for the test (default: 0.05).
#' @param num_simulations The number of simulations to run (default: 1000). More
#'   simulations increase accuracy but take longer.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" or "less".
#'
#' @return An estimated power value (proportion of simulations correctly rejecting H0).
#' @examples
#' # Example Data
#' set.seed(123) # for reproducibility
#' g1_data <- rnorm(30, mean = 5, sd = 1.5)
#' g2_data <- rnorm(40, mean = 6, sd = 2.5) # Different n and sd
#'
#' # Calculate power to detect a true difference of 1.5
#' power_est <- power_welch_simulation(g1_data, g2_data, delta = 1.5)
#' print(paste("Estimated Power (delta=1.5):", round(power_est, 3)))
#'
#' # Calculate power to detect a smaller difference of 0.5
#' power_est_small <- power_welch_simulation(g1_data, g2_data, delta = 0.5)
#' print(paste("Estimated Power (delta=0.5):", round(power_est_small, 3)))
#'
#' # One-sided test example (power to detect group 2 mean > group 1 mean by 1.0)
#' power_one_sided <- power_welch_simulation(g1_data, g2_data, delta = 1.0, alternative = "greater")
#' print(paste("Estimated Power (delta=1.0, one-sided):", round(power_one_sided, 3)))

power_welch_simulation <- function(group1_values,
                                   group2_values,
                                   delta,
                                   sig.level = 0.05,
                                   num_simulations = 1000,
                                   alternative = "two.sided") {

  # --- Input Validation ---
  if (!is.numeric(group1_values) || !is.numeric(group2_values)) {
    stop("group1_values and group2_values must be numeric vectors.")
  }
  if (length(group1_values) < 2 || length(group2_values) < 2) {
    stop("Both groups must have at least 2 observations to calculate variance.")
  }
  if (!is.numeric(delta) || length(delta) != 1) {
    stop("delta must be a single numeric value.")
  }
  if (!is.numeric(sig.level) || sig.level <= 0 || sig.level >= 1) {
    stop("sig.level must be between 0 and 1.")
  }
  if (!is.numeric(num_simulations) || num_simulations <= 0 || num_simulations != round(num_simulations)) {
    stop("num_simulations must be a positive integer.")
  }
  if (!alternative %in% c("two.sided", "less", "greater")) {
        stop("alternative must be one of 'two.sided', 'less', or 'greater'")
  }


  # --- Extract Characteristics from Data ---
  n1 <- length(group1_values)
  n2 <- length(group2_values)
  sd1 <- sd(group1_values)
  sd2 <- sd(group2_values)

  if (is.na(sd1) || is.na(sd2) || sd1 <= 0 || sd2 <= 0) {
      stop("Could not calculate valid positive standard deviations for both groups.")
  }

  # --- Simulation Loop ---
  message(paste("Running", num_simulations, "simulations..."))
  p_values <- numeric(num_simulations) # Store p-values from each simulation

  # Define means under H1 - we can set group 1 mean to 0 WLOG
  mean1_h1 <- 0
  mean2_h1 <- delta # Group 2 mean is delta away from group 1 mean

  for (i in 1:num_simulations) {
    # Simulate data under the alternative hypothesis (H1)
    sim_g1 <- rnorm(n1, mean = mean1_h1, sd = sd1)
    sim_g2 <- rnorm(n2, mean = mean2_h1, sd = sd2)

    # Perform Welch's t-test on the simulated data
    # Use suppressWarnings for potential minor issues in t.test under simulation
    test_result <- suppressWarnings(
        t.test(sim_g1, sim_g2, var.equal = FALSE, alternative = alternative)
    )

    p_values[i] <- test_result$p.value
  }

  # --- Calculate Power ---
  # Power is the proportion of simulations where the null hypothesis (H0) was rejected
  power_estimate <- mean(p_values < sig.level, na.rm = TRUE) # na.rm just in case a test fails catastrophically

  message("Simulation complete.")
  return(power_estimate)
}


simulation_results <-read.csv("DGP1_500_results.csv", row.names = 1)

sink("DGP1_500_power.txt")

mvbcf_cols <- c('mvbcf_tau_951','mvbcf_tau_952','mvbcf_tau_951','mvbcf_tau_952','mvbcf_pehe1','mvbcf_pehe2')

bcf_cols <- c('wsbcf_tau_951','wsbcf_tau_952','mvbart_tau_951','mvbart_tau_952','bcf_pehe1','bcf_pehe2')

for (j in 1:length(mvbcf_cols)) {

cat(paste(mvbcf_cols[j],"vs", bcf_cols[j], "\n"))

sample_group1 <- simulation_results[[mvbcf_cols[j]]]
sample_group2 <- simulation_results[[bcf_cols[j]]]
  
current_delta <- 0.015

# Check the condition for j
if (j == 5 || j == 6) {
  current_delta <- 0.5
}

power_val <- power_welch_simulation(
  group1_values = sample_group1,
  group2_values = sample_group2,
  delta = current_delta,        
  sig.level = 0.05,
  num_simulations = 500
)
print(paste("Estimated Power (delta=",current_delta,")",":", round(power_val, 3)))
cat("\n")

}



sink()