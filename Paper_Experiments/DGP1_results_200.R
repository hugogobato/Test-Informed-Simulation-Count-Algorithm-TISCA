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



#' Compare Multiple Metrics Between Two Models Using Robust Location Tests
#'
#' Iterates through pairs of corresponding metrics for two models, performs a robust
#' location test (user choice: Welch or Fligner-Policello) for each pair, and
#' applies multiple comparison correction across all location tests performed.
#' Scale comparison is NOT performed.
#'
#' @param data A data frame or list containing the simulation results.
#' @param model1_cols Character vector of column names in `data` for the first model's metrics.
#' @param model2_cols Character vector of column names in `data` for the second model's metrics.
#'        Must be the same length as `model1_cols` and correspond element-wise.
#' @param location_test_method Character string specifying the location test.
#'        Allowed values: "welch" (default) or "fligner-policello".
#' @param alpha Numeric significance level (default: 0.05).
#' @param correction_methods Character vector of p-value adjustment methods
#'        to apply (passed to `p.adjust.methods`). Default includes common ones.
#'
#' @return A list containing the results:
#'         - `parameters`: Input parameters used for the analysis.
#'         - `results_by_metric`: A list where each element corresponds to a metric pair,
#'           containing the raw location test object and original p-value for that metric.
#'         - `all_location_p_values_original`: Named vector of all original location p-values
#'           collected across metrics.
#'         - `p_values_adjusted`: A matrix of adjusted p-values (rows = location tests, cols = methods).
#'         - `significant_results`: A matrix/data frame indicating significance (TRUE/FALSE) based on adjusted p-values and alpha.
#'         - `warnings_list`: A list of warnings generated during execution for specific metrics.
#'
#' @examples
#' \dontrun{
#' # --- Example Data Setup ---
#' set.seed(123)
#' n_sim <- 50
#' simulation_results <- data.frame(
#'   mvbcf_pehe1 = rnorm(n_sim, 1.0, 0.2),
#'   mvbcf_ate1 = rnorm(n_sim, 0.5, 0.1),
#'   mvbcf_cover1 = rnorm(n_sim, 0.95, 0.05),
#'   bcf_pehe1 = rnorm(n_sim, 1.2, 0.25),  # Worse PEHE
#'   bcf_ate1 = rnorm(n_sim, 0.5, 0.12),   # Similar ATE, slightly more variance
#'   bcf_cover1 = rnorm(n_sim, 0.90, 0.05),  # Worse coverage
#'   mvbcf_runtime = rexp(n_sim, rate = 1/20), # seconds
#'   bcf_runtime = rexp(n_sim, rate = 1/15)   # Faster runtime
#' )
#' # Add some NAs
#' simulation_results$mvbcf_pehe1[sample(n_sim, 3)] <- NA
#' simulation_results$bcf_runtime[sample(n_sim, 2)] <- NA
#'
#' # Define corresponding columns
#' model1_cols <- c("mvbcf_pehe1", "mvbcf_ate1", "mvbcf_cover1", "mvbcf_runtime")
#' model2_cols <- c("bcf_pehe1", "bcf_ate1", "bcf_cover1", "bcf_runtime")
#'
#' # --- Run the Comparison ---
#'
#' # Using Welch's t-test
#' comparison_results_welch <- compare_multiple_locations(
#'   data = simulation_results,
#'   model1_cols = model1_cols,
#'   model2_cols = model2_cols,
#'   location_test_method = "welch"
#' )
#'
#' print(comparison_results_welch$p_values_adjusted)
#' print(comparison_results_welch$significant_results)
#' print(comparison_results_welch$warnings_list)
#'
#' # Using Fligner-Policello test
#' comparison_results_fp <- compare_multiple_locations(
#'   data = simulation_results,
#'   model1_cols = model1_cols,
#'   model2_cols = model2_cols,
#'   location_test_method = "fligner-policello"
#' )
#'
#' print(comparison_results_fp$p_values_adjusted)
#' print(comparison_results_fp$warnings_list)
#' }
#' 


compare_multiple_locations <- function(data,
                                       model1_cols,
                                       model2_cols,
                                       location_test_method = "welch",
                                       alpha = 0.05,
                                       correction_methods = c("bonferroni", "holm", "BH")) {

  # --- Input Validation ---
  if (!is.data.frame(data) && !is.list(data)) {
    stop("'data' must be a data frame or a list.", call. = FALSE)
  }
  if (length(model1_cols) != length(model2_cols)) {
    stop("'model1_cols' and 'model2_cols' must have the same length.", call. = FALSE)
  }
  if (length(model1_cols) == 0) {
    stop("Column lists cannot be empty.", call. = FALSE)
  }
  if (!all(sapply(model1_cols, function(cn) cn %in% names(data)))) {
      missing_cols <- model1_cols[!model1_cols %in% names(data)]
      stop("Missing columns for model 1 in data: ", paste(missing_cols, collapse=", "), call. = FALSE)
  }
  if (!all(sapply(model2_cols, function(cn) cn %in% names(data)))) {
      missing_cols <- model2_cols[!model2_cols %in% names(data)]
      stop("Missing columns for model 2 in data: ", paste(missing_cols, collapse=", "), call. = FALSE)
  }
  if (!location_test_method %in% c("welch", "fligner-policello")) {
    stop("location_test_method must be 'welch' or 'fligner-policello'.", call. = FALSE)
  }
  if (location_test_method == "fligner-policello" && !requireNamespace("kSamples", quietly = TRUE)) {
     stop("Package 'kSamples' needed for Fligner-Policello test. Please install it.", call. = FALSE)
  }
  # No longer need 'car' package
  # if (!requireNamespace("car", quietly = TRUE)) {
  #  stop("Package 'car' needed for Levene's test. Please install it.", call. = FALSE)
  #}
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1.", call. = FALSE)
  }

  # --- Initialization ---
  all_location_p_values_list <- list() # Collect all location p-values
  results_by_metric <- list()          # Store detailed results per metric
  warnings_list <- list()              # Collect warnings

  num_metrics <- length(model1_cols)
  message(paste("Starting location comparison for", num_metrics, "metric(s)..."))

  # --- Loop Through Metrics ---
  for (i in 1:num_metrics) {
    col1 <- model1_cols[i]
    col2 <- model2_cols[i]
    # Generate a more robust metric name base
    split1 <- strsplit(col1, "_")[[1]]
    split2 <- strsplit(col2, "_")[[1]]
    prefix_len <- 0
    if (length(split1) > 1 && length(split2) > 1 && split1[1] == split2[1]){
        prefix_len <- nchar(split1[1]) + 1 # length of prefix + underscore
    }
    metric_name <- substr(col1, prefix_len + 1, nchar(col1))
    if (nchar(metric_name) == 0) metric_name <- paste0("metric_", i) # Fallback

    message(paste("Processing metric:", metric_name, "(", col1, "vs", col2, ")"))

    metric_results <- list(metric_name = metric_name, model1_col = col1, model2_col = col2)
    metric_warnings <- list()

    # Extract and clean data for the current metric
    group1_raw <- data[[col1]]
    group2_raw <- data[[col2]]

    if(is.list(group1_raw)) group1_raw <- unlist(group1_raw)
    if(is.list(group2_raw)) group2_raw <- unlist(group2_raw)

    group1 <- tryCatch(as.numeric(group1_raw), warning = function(w) {
        metric_warnings <<- c(metric_warnings, paste("Coercion warning for", col1, ":", w$message))
        suppressWarnings(as.numeric(group1_raw))
    })
    group2 <- tryCatch(as.numeric(group2_raw), warning = function(w) {
        metric_warnings <<- c(metric_warnings, paste("Coercion warning for", col2, ":", w$message))
        suppressWarnings(as.numeric(group2_raw))
    })

    # Handle NAs
    valid_idx1 <- !is.na(group1)
    valid_idx2 <- !is.na(group2)
    group1_complete <- group1[valid_idx1]
    group2_complete <- group2[valid_idx2]

    n1 <- length(group1_complete)
    n2 <- length(group2_complete)

    if (sum(!valid_idx1) > 0) metric_warnings <- c(metric_warnings, paste(sum(!valid_idx1), "NA(s) removed from", col1))
    if (sum(!valid_idx2) > 0) metric_warnings <- c(metric_warnings, paste(sum(!valid_idx2), "NA(s) removed from", col2))

    # Initialize p-value for this metric as NA
    p_loc <- NA
    loc_test_name_suffix <- ""
    location_test_result <- NULL

    # --- Perform Location Test if data sufficient ---
    if (n1 < 2 || n2 < 2) { # t.test needs at least 2, FP often needs more but min check is 2
      metric_warnings <- c(metric_warnings, "Insufficient non-NA data (< 2 in at least one group) for location testing.")
      message("  Skipping location test due to insufficient data.")
    } else {
      # Location Test
      loc_test_failed <- FALSE
      if (location_test_method == "welch") {
        loc_test_name_suffix <- "_welch"
        location_test_result <- tryCatch({
          t.test(group1_complete, group2_complete, var.equal = FALSE)
        }, error = function(e) {
          metric_warnings <<- c(metric_warnings, paste("Welch t-test failed:", e$message))
          loc_test_failed <<- TRUE
          NULL
        })
        if(!is.null(location_test_result)) p_loc <- location_test_result$p.value

      } else if (location_test_method == "fligner-policello") {
        loc_test_name_suffix <- "_fp"
        min_n_fp <- 10 # Recommended minimum sample size per group
        if (n1 < min_n_fp || n2 < min_n_fp) {
             warn_msg <- paste0("Sample size(s) (", n1, ", ", n2,
                             ") may be small for reliable Fligner-Policello results (recommend >=", min_n_fp, ").")
             metric_warnings <- c(metric_warnings, warn_msg)
             message(paste("  Warning:", warn_msg))
         }
         location_test_result <- tryCatch({
            # Check if data is constant, fp.test might fail
            if (length(unique(group1_complete)) <= 1 || length(unique(group2_complete)) <= 1) {
                stop("Data within at least one group is constant.")
            }
            fp.test(group1_complete, group2_complete)
         }, error = function(e) {
            metric_warnings <<- c(metric_warnings, paste("Fligner-Policello test failed:", e$message))
            loc_test_failed <<- TRUE
            NULL
         })
        if(!is.null(location_test_result)) p_loc <- location_test_result$p.value
      }
       metric_results$location_test_result <- location_test_result

    } # End if sufficient data

    # Store location p-value with unique name
    p_loc_name <- paste0(metric_name, "_location", loc_test_name_suffix)
    all_location_p_values_list[[p_loc_name]] <- p_loc

    # Store metric-specific results
    metric_results$p_value_location_original <- p_loc
    metric_results$warnings <- metric_warnings
    results_by_metric[[metric_name]] <- metric_results
    if(length(metric_warnings) > 0) warnings_list[[metric_name]] <- metric_warnings

  } # End loop through metrics

  message("Finished individual location tests. Applying multiple comparison corrections...")

  # --- Apply Multiple Comparison Correction ---
  p_values_vec <- unlist(all_location_p_values_list)
  valid_p_indices <- !is.na(p_values_vec)
  p_values_to_adjust <- p_values_vec[valid_p_indices]

  p_adj_matrix <- NULL
  sig_matrix <- NULL

  if (length(p_values_to_adjust) > 0) {
      # Ensure only valid methods are requested
      valid_methods <- p.adjust.methods
      correction_methods <- intersect(correction_methods, valid_methods)
       if(length(correction_methods) == 0) {
          warning("No valid correction methods specified or available. Using 'bonferroni'.", call.=FALSE)
          correction_methods <- "bonferroni"
      }

      p_adj_list <- lapply(correction_methods, function(method) {
          p.adjust(p_values_to_adjust, method = method)
      })
      p_adj_matrix_valid <- do.call(cbind, p_adj_list)
      colnames(p_adj_matrix_valid) <- correction_methods
      rownames(p_adj_matrix_valid) <- names(p_values_to_adjust)

      # Create full matrix including NAs for alignment
      p_adj_matrix <- matrix(NA_real_, nrow = length(p_values_vec), ncol = length(correction_methods),
                            dimnames = list(names(p_values_vec), correction_methods))
      p_adj_matrix[valid_p_indices, ] <- p_adj_matrix_valid

      # Determine significance based on adjusted values
      sig_matrix <- p_adj_matrix < alpha
      # Ensure NAs remain NA in significance matrix
      sig_matrix[is.na(p_adj_matrix)] <- NA

      message(paste("Corrections applied to", length(p_values_to_adjust), "valid p-values."))

  } else {
      message("No valid p-values obtained across all metrics to apply corrections.")
  }


  # --- Consolidate Results ---
  final_results <- list(
      parameters = list(
          location_test_method = location_test_method,
          alpha = alpha,
          correction_methods = correction_methods,
          num_metrics = num_metrics,
          num_tests_corrected = length(p_values_to_adjust)
      ),
      results_by_metric = results_by_metric,
      all_location_p_values_original = p_values_vec, # Vector with NAs included
      p_values_adjusted = p_adj_matrix,
      significant_results = sig_matrix,
      warnings_list = warnings_list
  )

  message("Comparison complete.")
  return(final_results)
}

simulation_results <-read.csv("DGP1_200_results.csv", row.names = 1)

sink("DGP1_200.txt")

mvbcf_cols <- c('mvbcf_tau_951','mvbcf_tau_952','mvbcf_pehe1','mvbcf_pehe2','mvbart_tau_951','mvbart_tau_952')

bcf_cols <- c('wsbcf_tau_951','wsbcf_tau_952','bcf_pehe1','bcf_pehe2','mvbcf_tau_951','mvbcf_tau_952')

cat("\n\n--- Running Example with Welch (Location Only) ---\n")
comparison_results_welch <- compare_multiple_locations(
   data = simulation_results,
   model1_cols =mvbcf_cols,
   model2_cols = bcf_cols,
   location_test_method = "welch"
 )

cat("\n--- Original P-values (Welch) ---\n")
print(comparison_results_welch$all_location_p_values_original)
cat("\n--- Adjusted P-values (Welch) ---\n")
print(comparison_results_welch$p_values_adjusted)
cat("\n--- Significance (Welch, alpha=0.05) ---\n")
print(comparison_results_welch$significant_results)
cat("\n--- Warnings (Welch) ---\n")
print(comparison_results_welch$warnings_list)


sink()