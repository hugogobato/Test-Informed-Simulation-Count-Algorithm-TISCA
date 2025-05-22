################################################################################
# TISCA (Test-Informed Simulation Count Algorithm) Implementation
# Integrating simulation code with TISCA methodology
################################################################################

# Load required libraries
library(progress)
library(mvtnorm)
library(stochtree)
library(mvbcf)
library(dbarts)
library(bartCause)
library(skewBART)
library(ggplot2)
library(reshape2)

################################################################################
# Helper Functions
################################################################################

# Function to check if a value is within a credible interval
in_cred <- function(samples, value, interval) {
  upper_quantile <- 1-(1-interval)/2
  lower_quantile <- 0+(1-interval)/2
  
  q1 <- quantile(samples, lower_quantile)
  q2 <- quantile(samples, upper_quantile)
  
  in_cred <- ifelse(value >= q1 & value <= q2, T, F)
}

# Function to calculate credible interval width
cred_width <- function(samples, interval) {
  upper_quantile <- 1-(1-interval)/2
  lower_quantile <- 0+(1-interval)/2
  
  q1 <- quantile(samples, lower_quantile)
  q2 <- quantile(samples, upper_quantile)
  
  return(q2-q1)
}

# Helper function for Welch's t-test
welch_t_test <- function(data, metric_p, metric_b) {
  # Extract the metrics for both models
  values_p <- data[[metric_p]]
  values_b <- data[[metric_b]]
  
  # Perform Welch's t-test
  test_result <- t.test(values_p, values_b, paired = FALSE, var.equal = FALSE)
  
  return(list(
    p_value = test_result$p.value,
    t_stat = test_result$statistic,
    sd_p = sd(values_p),
    sd_b = sd(values_b),
    mean_p = mean(values_p),
    mean_b = mean(values_b)
  ))
}

# Helper function to adjust p-values for multiple testing
adjust_p_values <- function(p_values, method) {
  if (method == "none") {
    return(p_values)
  } else {
    return(p.adjust(p_values, method = method))
  }
}

# Helper function to estimate power for Welch's t-test
estimate_welch_power <- function(n1, n2, sd1, sd2, delta, alpha) {
  # Calculate degrees of freedom using Welch-Satterthwaite equation
  df <- ((sd1^2/n1 + sd2^2/n2)^2) / 
        ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
  
  # Calculate non-centrality parameter
  ncp <- delta / sqrt(sd1^2/n1 + sd2^2/n2)
  
  # Calculate critical value
  crit <- qt(1-alpha/2, df)
  
  # Calculate power
  power <- 1 - pt(crit, df, ncp) + pt(-crit, df, ncp)
  
  return(power)
}

################################################################################
# Simulation Function
################################################################################

# Define the simulation function that will be called by TISCA
run_simulation <- function(seed) {
  # Set random seed
  set.seed(seed)
  
  # Train Data
  n <- 500
  
  X1 <- runif(n)
  X2 <- runif(n)
  X3 <- runif(n)
  X4 <- runif(n)
  X5 <- runif(n)
  X6 <- rbinom(n, 1, 0.5)
  X7 <- rbinom(n, 1, 0.5)
  X8 <- rbinom(n, 1, 0.5)
  X9 <- sample(c(0, 1, 2, 3, 4), n, replace=T)
  X10 <- sample(c(0, 1, 2, 3, 4), n, replace=T)
  
  X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  
  Mu1 <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
  Mu2 <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
  
  Tau1 <- (2*X4+2*X5)*10
  Tau2 <- (1*X4+3*X5)*10
  
  true_propensity <- X4
  
  Z <- rbinom(n, 1, true_propensity)
  
  Y <- cbind(Mu1+Z*Tau1, Mu2+Z*Tau2) + mvtnorm::rmvnorm(n, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow=2, byrow=T))
  
  # Test Data
  n_test <- 1000
  
  X1_test <- runif(n_test)
  X2_test <- runif(n_test)
  X3_test <- runif(n_test)
  X4_test <- runif(n_test)
  X5_test <- runif(n_test)
  X6_test <- rbinom(n_test, 1, 0.5)
  X7_test <- rbinom(n_test, 1, 0.5)
  X8_test <- rbinom(n_test, 1, 0.5)
  X9_test <- sample(c(0, 1, 2, 3, 4), n_test, replace=T)
  X10_test <- sample(c(0, 1, 2, 3, 4), n_test, replace=T)
  
  X_test <- cbind(X1_test, X2_test, X3_test, X4_test, X5_test, X6_test, X7_test, X8_test, X9_test, X10_test)
  
  Mu1_test <- (11*sin(pi*X1_test*X2_test)+18*(X3_test-0.5)^2+10*X4_test+12*X6_test+X9_test)*10+300
  Mu2_test <- (9*sin(pi*X1_test*X2_test)+22*(X3_test-0.5)^2+14*X4_test+8*X6_test+X9_test)*10+300
  
  Tau1_test <- (2*X4_test+2*X5_test)*10
  Tau2_test <- (1*X4_test+3*X5_test)*10
  
  true_propensity_test <- X4_test
  
  Z_test <- rbinom(n_test, 1, true_propensity_test)
  
  Y_test <- cbind(Mu1_test+Z_test*Tau1_test, Mu2_test+Z_test*Tau2_test) + mvtnorm::rmvnorm(n_test, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow=2, byrow=T))
  
  # Estimate of propensity score
  p_mod <- bart(x.train = X, y.train = Z, x.test = X_test, k=3, verbose = FALSE)
  p <- colMeans(pnorm(p_mod$yhat.train))
  p_test <- colMeans(pnorm(p_mod$yhat.test))
  
  # Adding to matrix
  X2 <- X
  X2_test <- X_test
  X <- cbind(X, p)
  X_test <- cbind(X_test, p_test)
  Z2 <- cbind(Z,Z)
  
  # Set some parameters
  n_tree_mu <- 50
  n_tree_tau <- 20
  n_iter <- 1000
  n_burn <- 500
  
  mu_val <- 1
  tau_val <- 0.375
  v_val <- 1
  wish_val <- 1
  min_val <- 1
  
  # MVBCF model
  mvbcf_mod <- run_mvbcf(X,
                         Y,
                         Z,
                         X2,
                         X_test, # Test data for the mu part of the model
                         X2_test, # Test data for the tau part of the model
                         0.95,
                         2,
                         0.25,
                         3,
                         diag((mu_val)^2/n_tree_mu, 2),
                         diag((tau_val)^2/n_tree_tau, 2),
                         v_val,
                         diag(wish_val, 2),
                         n_iter,
                         n_tree_mu,
                         n_tree_tau,
                         min_val)
  
  mvbcf_tau_preds1 <- rowMeans(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)])
  mvbcf_ate1 <- mean(mvbcf_tau_preds1)
  mvbcf_tau_preds2 <- rowMeans(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)])
  mvbcf_ate2 <- mean(mvbcf_tau_preds2)
  
  mvbcf_pehe1 <- sqrt(mean((Tau1_test-mvbcf_tau_preds1)^2))
  mvbcf_pehe2 <- sqrt(mean((Tau2_test-mvbcf_tau_preds2)^2))
  
  # BCF model 1
  bcfmod1 <- bcf(X2, Z, Y[,1], X_test = X2_test,
                prognostic_forest_params=list("num_trees"= n_tree_mu, "sigma2_leaf_init"=sd(Y[,1])^2),
                treatment_effect_forest_params=list("num_trees"= n_tree_tau, "sigma2_leaf_init"=(0.375*sd(Y[,1]))^2),
                Z_test = Z_test, propensity_test = p_test,
                propensity_train=p, num_burnin=n_burn, num_mcmc=n_iter-n_burn, num_gfr = 0)
  bcf_pehe1 <- sqrt(mean((Tau1_test-rowMeans(bcfmod1$tau_hat_test))^2))
  bcf_ate1 <- mean(rowMeans(bcfmod1$tau_hat_test))
  
  # BCF model 2
  bcfmod2 <- bcf(X2, Z, Y[,2], X_test = X2_test,
                prognostic_forest_params=list("num_trees"= n_tree_mu, "sigma2_leaf_init"=sd(Y[,1])^2),
                treatment_effect_forest_params=list("num_trees"= n_tree_tau, "sigma2_leaf_init"=(0.375*sd(Y[,1]))^2),
                Z_test = Z_test, propensity_test = p_test,
                propensity_train=p, num_burnin=n_burn, num_mcmc=n_iter-n_burn, num_gfr = 0)
  bcf_pehe2 <- sqrt(mean((Tau2_test-rowMeans(bcfmod2$tau_hat_test))^2))
  bcf_ate2 <- mean(rowMeans(bcfmod2$tau_hat_test))
  
  # WSBCF model 1
  wsbcfmod1 <- bcf(X2, Z, Y[,1], X_test = X2_test,
                  Z_test = Z_test, propensity_test = p_test,
                  propensity_train=p, num_burnin=n_burn, num_mcmc=n_iter-n_burn, num_gfr = 50)
  wsbcf_pehe1 <- sqrt(mean((Tau1_test-rowMeans(wsbcfmod1$tau_hat_test))^2))
  wsbcf_ate1 <- mean(rowMeans(wsbcfmod1$tau_hat_test))
  
  # WSBCF model 2
  wsbcfmod2 <- bcf(X2, Z, Y[,2], X_test = X2_test,
                  Z_test = Z_test, propensity_test = p_test,
                  propensity_train=p, num_burnin=n_burn, num_mcmc=n_iter-n_burn, num_gfr = 50)
  wsbcf_pehe2 <- sqrt(mean((Tau2_test-rowMeans(wsbcfmod2$tau_hat_test))^2))
  wsbcf_ate2 <- mean(rowMeans(bcfmod2$tau_hat_test))
  
  ate1 <- mean(Tau1_test)
  ate2 <- mean(Tau2_test)
  
  # Calculate coverage metrics
  mvbcf_tau_951 <- mean(diag(apply(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)], 1, in_cred, Tau1_test, 0.95)))
  mvbcf_tau_951w <- mean(apply(mvbcf_mod$predictions_tau_test[,1,-c(1:n_burn)], 1, cred_width, 0.95))
  bcf_tau_951 <- mean(diag(apply(bcfmod1$tau_hat_test, 1, in_cred, Tau1_test, 0.95)))
  bcf_tau_951w <- mean(apply(bcfmod1$tau_hat_test, 1, cred_width, 0.95))
  
  mvbcf_tau_952 <- mean(diag(apply(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)], 1, in_cred, Tau2_test, 0.95)))
  mvbcf_tau_952w <- mean(apply(mvbcf_mod$predictions_tau_test[,2,-c(1:n_burn)], 1, cred_width, 0.95))
  bcf_tau_952 <- mean(diag(apply(bcfmod2$tau_hat_test, 1, in_cred, Tau2_test, 0.95)))
  bcf_tau_952w <- mean(apply(bcfmod2$tau_hat_test, 1, cred_width, 0.95))
  
  wsbcf_tau_951 <- mean(diag(apply(wsbcfmod1$tau_hat_test, 1, in_cred, Tau1_test, 0.95)))
  wsbcf_tau_952 <- mean(diag(apply(wsbcfmod2$tau_hat_test, 1, in_cred, Tau2_test, 0.95)))
  
  # BART models
  colnames(X) <- rep("", ncol(X))
  colnames(X_test) <- rep("", ncol(X_test))
  colnames(X) <- paste0("V", 1:ncol(X))
  colnames(X_test) <- paste0("V", 1:ncol(X_test))
  
  bart_mod1 <- bartc(Y[,1], Z, X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau)
  bart_tau_preds1 <- colMeans(predict(bart_mod1, X_test, type="icate"))
  bart_ate1 <- mean(bart_tau_preds1)
  
  bart_mod2 <- bartc(Y[,2], Z, X, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_tau)
  bart_tau_preds2 <- colMeans(predict(bart_mod2, X_test, type="icate"))
  bart_ate2 <- mean(bart_tau_preds2)
  
  # Multivariate skew BART
  hypers <- Hypers(X = cbind(X, Z), Y = Y, num_tree = n_tree_mu+n_tree_tau)
  opts <- Opts(num_burn = n_burn, num_save = n_iter-n_burn)
  X_test_z <- rbind(cbind(X_test, 0*Z_test), cbind(X_test, 0*Z_test+1))
  
  fitted_Multiskewbart <- MultiskewBART(X = cbind(X, Z), Y = Y, test_X = X_test_z, hypers=hypers, opts=opts, do_skew = F)
  
  z_0_preds <- fitted_Multiskewbart$y_hat_test[1:n_test,,]
  z_1_preds <- fitted_Multiskewbart$y_hat_test[-c(1:n_test),,]
  z_1_0_preds <- z_1_preds-z_0_preds
  
  mvbart_icates1 <- rowMeans(z_1_0_preds[,1,])
  mvbart_ate1 <- mean(mvbart_icates1)
  mvbart_icates2 <- rowMeans(z_1_0_preds[,2,])
  mvbart_ate2 <- mean(mvbart_icates2)
  
  bart_pehe1 <- sqrt(mean((Tau1_test-bart_tau_preds1)^2))
  mvbart_pehe1 <- sqrt(mean((Tau1_test-mvbart_icates1)^2))
  bart_tau_951 <- mean(diag(apply(predict(bart_mod1, X_test, type="icate"), 2, in_cred, Tau1_test, 0.95)))
  mvbart_tau_951 <- mean(diag(apply(z_1_0_preds[,1,], 1, in_cred, Tau1_test, 0.95)))
  
  bart_pehe2 <- sqrt(mean((Tau2_test-bart_tau_preds2)^2))
  mvbart_pehe2 <- sqrt(mean((Tau2_test-mvbart_icates2)^2))
  bart_tau_952 <- mean(diag(apply(predict(bart_mod2, X_test, type="icate"), 2, in_cred, Tau2_test, 0.95)))
  mvbart_tau_952 <- mean(diag(apply(z_1_0_preds[,2,], 1, in_cred, Tau2_test, 0.95)))
  
  # Return results as a data frame
  results <- data.frame(
    mvbcf_pehe1 = mvbcf_pehe1,
    mvbcf_pehe2 = mvbcf_pehe2,
    bcf_pehe1 = bcf_pehe1,
    bcf_pehe2 = bcf_pehe2,
    wsbcf_pehe1 = wsbcf_pehe1,
    wsbcf_pehe2 = wsbcf_pehe2,
    bart_pehe1 = bart_pehe1,
    bart_pehe2 = bart_pehe2,
    mvbart_pehe1 = mvbart_pehe1,
    mvbart_pehe2 = mvbart_pehe2,
    mvbcf_tau_951 = mvbcf_tau_951,
    mvbcf_tau_952 = mvbcf_tau_952,
    bcf_tau_951 = bcf_tau_951,
    bcf_tau_952 = bcf_tau_952,
    wsbcf_tau_951 = wsbcf_tau_951,
    wsbcf_tau_952 = wsbcf_tau_952,
    bart_tau_951 = bart_tau_951,
    bart_tau_952 = bart_tau_952,
    mvbart_tau_951 = mvbart_tau_951,
    mvbart_tau_952 = mvbart_tau_952
  )
  
  return(results)
}

################################################################################
# TISCA Main Function
################################################################################

run_tisca <- function(
  sim_func,                # Function that runs one simulation and returns metrics
  comparison_pairs,        # List of pairs to compare (metric_p, metric_b)
  mdes,                    # Vector of minimum detectable effect sizes
  target_power = 0.80,     # Target statistical power
  alpha = 0.05,            # Significance level
  batch_size = 50,         # Number of simulations per batch
  initial_count = NULL,    # Initial simulation count (defaults to batch_size if NULL)
  correction_method = "none", # Multiple testing correction method
  verbose = TRUE,          # Whether to print progress updates
  save_results = TRUE      # Whether to save results to CSV files
) {
  # Set initial count to batch_size if not specified
  if (is.null(initial_count)) {
    initial_count <- batch_size
  }
  
  # Initialization Phase
  J <- 0
  results_agg <- data.frame() # Empty dataframe to store results
  K <- length(comparison_pairs)
  P_current <- rep(0, K)      # Vector of K zeros for current power
  
  # Create data structures to track power and p-values over iterations
  power_tracking <- data.frame(
    simulations = numeric(),
    comparison = character(),
    power = numeric()
  )
  
  pvalue_tracking <- data.frame(
    simulations = numeric(),
    comparison = character(),
    p_raw = numeric(),
    p_adj = numeric()
  )
  
  # Create a progress bar if verbose
  if (verbose) {
    cat("Starting TISCA algorithm...\n")
    cat("Target power:", target_power, "\n")
    cat("Significance level:", alpha, "\n")
    cat("Batch size:", batch_size, "\n")
    cat("Initial count:", initial_count, "\n")
    cat("Correction method:", correction_method, "\n")
    cat("Number of comparisons:", K, "\n\n")
  }
  
  # Optional Initial Run
  if (initial_count > 0 && initial_count >= batch_size) {
    if (verbose) {
      cat("Running initial", initial_count, "simulations...\n")
      pb <- progress_bar$new(total = initial_count)
    }
    
    # Run initial simulations
    for (i in 1:initial_count) {
      if (verbose) pb$tick()
      metrics_run <- sim_func(seed = i)
      results_agg <- rbind(results_agg, metrics_run)
    }
    J <- initial_count
    
    # Perform initial power calculation
    p_values_raw <- numeric(K)
    test_stats <- numeric(K)
    mean_diffs <- numeric(K)
    
    for (k in 1:K) {
      # Extract comparison information
      metric_p <- comparison_pairs[[k]][1]
      metric_b <- comparison_pairs[[k]][2]
      comparison_name <- paste(metric_p, "vs", metric_b)
      
      # Perform Welch's t-test
      test_result <- welch_t_test(results_agg, metric_p, metric_b)
      
      # Store results
      p_values_raw[k] <- test_result$p_value
      test_stats[k] <- test_result$t_stat
      mean_diffs[k] <- test_result$mean_p - test_result$mean_b
      
      # Estimate power if standard deviations are valid
      if (test_result$sd_p > 0 && test_result$sd_b > 0) {
        P_current[k] <- estimate_welch_power(
          J, J, 
          test_result$sd_p, test_result$sd_b, 
          mdes[k], alpha
        )
      } else {
        P_current[k] <- 0
      }
      
      # Track power
      power_tracking <- rbind(power_tracking, data.frame(
        simulations = J,
        comparison = comparison_name,
        power = P_current[k]
      ))
    }
    
    # Adjust p-values for multiple testing
    p_values_adj <- adjust_p_values(p_values_raw, correction_method)
    
    # Track p-values
    for (k in 1:K) {
      comparison_name <- paste(comparison_pairs[[k]][1], "vs", comparison_pairs[[k]][2])
      pvalue_tracking <- rbind(pvalue_tracking, data.frame(
        simulations = J,
        comparison = comparison_name,
        p_raw = p_values_raw[k],
        p_adj = p_values_adj[k]
      ))
    }
    
    if (verbose) {
      cat("\nAfter initial", J, "simulations:\n")
      for (k in 1:K) {
        cat(sprintf("Comparison %d: Power = %.4f, p-value = %.4f, adjusted p-value = %.4f\n", 
                    k, P_current[k], p_values_raw[k], p_values_adj[k]))
      }
      cat("Minimum power:", min(P_current), "\n\n")
    }
  }
  
  # Iterative Simulation and Power Check Loop
  iteration <- 1
  
  while (min(P_current) < target_power) {
    if (verbose) {
      cat(sprintf("Iteration %d: Running batch of %d simulations (total will be %d)...\n", 
                  iteration, batch_size, J + batch_size))
      pb <- progress_bar$new(total = batch_size)
    }
    
    # Run a new batch of simulations
    start_seed <- J + 1
    end_seed <- J + batch_size
    
    for (i in start_seed:end_seed) {
      if (verbose) pb$tick()
      metrics_run <- sim_func(seed = i)
      results_agg <- rbind(results_agg, metrics_run)
    }
    
    J <- J + batch_size
    
    # Perform tests and estimate power on current J simulations
    p_values_raw <- numeric(K)
    test_stats <- numeric(K)
    mean_diffs <- numeric(K)
    
    for (k in 1:K) {
      # Extract comparison information
      metric_p <- comparison_pairs[[k]][1]
      metric_b <- comparison_pairs[[k]][2]
      comparison_name <- paste(metric_p, "vs", metric_b)
      
      # Perform Welch's t-test
      test_result <- welch_t_test(results_agg, metric_p, metric_b)
      
      # Store results
      p_values_raw[k] <- test_result$p_value
      test_stats[k] <- test_result$t_stat
      mean_diffs[k] <- test_result$mean_p - test_result$mean_b
      
      # Estimate power if standard deviations are valid
      if (test_result$sd_p > 0 && test_result$sd_b > 0) {
        P_current[k] <- estimate_welch_power(
          J, J, 
          test_result$sd_p, test_result$sd_b, 
          mdes[k], alpha
        )
      } else {
        P_current[k] <- 0
      }
      
      # Track power
      power_tracking <- rbind(power_tracking, data.frame(
        simulations = J,
        comparison = comparison_name,
        power = P_current[k]
      ))
    }
    
    # Adjust p-values for multiple testing
    p_values_adj <- adjust_p_values(p_values_raw, correction_method)
    
    # Track p-values
    for (k in 1:K) {
      comparison_name <- paste(comparison_pairs[[k]][1], "vs", comparison_pairs[[k]][2])
      pvalue_tracking <- rbind(pvalue_tracking, data.frame(
        simulations = J,
        comparison = comparison_name,
        p_raw = p_values_raw[k],
        p_adj = p_values_adj[k]
      ))
    }
    
    if (verbose) {
      cat("\nAfter", J, "simulations:\n")
      for (k in 1:K) {
        cat(sprintf("Comparison %d: Power = %.4f, p-value = %.4f, adjusted p-value = %.4f\n", 
                    k, P_current[k], p_values_raw[k], p_values_adj[k]))
      }
      cat("Minimum power:", min(P_current), "\n\n")
    }
    
    iteration <- iteration + 1
  }
  
  # Output Phase
  J_final <- J
  P_raw <- p_values_raw
  P_adj <- p_values_adj
  Stats <- test_stats
  P_achieved <- P_current
  
  # Create summary table
  summary_table <- data.frame(
    Comparison = character(),
    Proposed_Metric = character(),
    Benchmark_Metric = character(),
    MDE = numeric(),
    Mean_Difference = numeric(),
    Power = numeric(),
    P_Value = numeric(),
    Adjusted_P_Value = numeric(),
    T_Statistic = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (k in 1:length(comparison_pairs)) {
    summary_table[k, "Comparison"] <- paste(comparison_pairs[[k]][1], "vs", comparison_pairs[[k]][2])
    summary_table[k, "Proposed_Metric"] <- comparison_pairs[[k]][1]
    summary_table[k, "Benchmark_Metric"] <- comparison_pairs[[k]][2]
    summary_table[k, "MDE"] <- mdes[k]
    summary_table[k, "Mean_Difference"] <- mean_diffs[k]
    summary_table[k, "Power"] <- P_achieved[k]
    summary_table[k, "P_Value"] <- P_raw[k]
    summary_table[k, "Adjusted_P_Value"] <- P_adj[k]
    summary_table[k, "T_Statistic"] <- Stats[k]
  }
  
  if (verbose) {
    cat("\n=== TISCA COMPLETED ===\n")
    cat("Final number of simulations required:", J_final, "\n")
    cat("Final achieved powers:\n")
    for (k in 1:K) {
      cat(sprintf("Comparison %d (%s vs %s): Power = %.4f, p-value = %.4f, adjusted p-value = %.4f\n", 
                  k, comparison_pairs[[k]][1], comparison_pairs[[k]][2],
                  P_achieved[k], P_raw[k], P_adj[k]))
    }
  }
  
  # Create power vs simulation plot
  if (save_results) {
    # Plot power vs simulations
    power_plot <- ggplot(power_tracking, aes(x = simulations, y = power, color = comparison, group = comparison)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      labs(
        title = "Power vs Number of Simulations",
        subtitle = paste("Target Power =", target_power),
        x = "Number of Simulations",
        y = "Estimated Power",
        color = "Comparison"
      ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    # Save power plot
    ggsave("power_vs_simulations.png", power_plot, width = 10, height = 6, dpi = 300)
    
    # Create p-value bar plot
    # Get the final p-values
    final_pvalues <- pvalue_tracking[pvalue_tracking$simulations == max(pvalue_tracking$simulations), ]
    
    # Transform p-values to -log10 scale for better visualization
    final_pvalues$neg_log10_p_raw <- -log10(final_pvalues$p_raw)
    final_pvalues$neg_log10_p_adj <- -log10(final_pvalues$p_adj)
    
    # Reshape for plotting
    pvalue_plot_data <- data.frame(
      comparison = rep(final_pvalues$comparison, 2),
      p_type = c(rep("Raw p-value", nrow(final_pvalues)), rep("Adjusted p-value", nrow(final_pvalues))),
      neg_log10_p = c(final_pvalues$neg_log10_p_raw, final_pvalues$neg_log10_p_adj)
    )
    
    # Create the bar plot
    pvalue_plot <- ggplot(pvalue_plot_data, aes(x = comparison, y = neg_log10_p, fill = p_type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      labs(
        title = "Significance of Comparisons",
        subtitle = paste("After", J_final, "Simulations"),
        x = "Comparison",
        y = "-log10(p-value)",
        fill = "P-value Type"
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      annotate("text", x = 1, y = -log10(0.05) + 0.2, label = "α = 0.05", color = "red")
    
    # Save p-value plot
    ggsave("pvalue_comparison.png", pvalue_plot, width = 10, height = 6, dpi = 300)
    
    # Save results to CSV
    write.csv(results_agg, "tisca_simulation_results.csv", row.names = FALSE)
    write.csv(summary_table, "tisca_summary_results.csv", row.names = FALSE)
    write.csv(power_tracking, "power_tracking.csv", row.names = FALSE)
    write.csv(pvalue_tracking, "pvalue_tracking.csv", row.names = FALSE)
    
    if (verbose) {
      cat("\nResults saved to CSV files\n")
      cat("Power vs Simulations plot saved as 'power_vs_simulations.png'\n")
      cat("P-value comparison plot saved as 'pvalue_comparison.png'\n")
    }
  }
  
  # Return results
  return(list(
    J_final = J_final,
    P_raw = P_raw,
    P_adj = P_adj,
    Stats = Stats,
    P_achieved = P_achieved,
    results = results_agg,
    comparison_pairs = comparison_pairs,
    mdes = mdes,
    summary_table = summary_table,
    power_tracking = power_tracking,
    pvalue_tracking = pvalue_tracking
  ))
}

################################################################################
# Example Usage
################################################################################

# Define comparison pairs (proposed model vs benchmark)
# Format: list(c("proposed_metric", "benchmark_metric"), ...)
comparison_pairs <- list(
  c("mvbcf_pehe1", "bcf_pehe1"),     # Compare MVBCF vs BCF for PEHE1
  c("mvbcf_pehe2", "bcf_pehe2"),     # Compare MVBCF vs BCF for PEHE2
  c("mvbcf_tau_951", "bcf_tau_951"), # Compare MVBCF vs BCF for coverage1
  c("mvbcf_tau_952", "bcf_tau_952")  # Compare MVBCF vs BCF for coverage2
)

# Define minimum detectable effect sizes
# Negative values for metrics where lower is better (like PEHE)
# Positive values for metrics where higher is better (like coverage)
mdes <- c(-0.5, -0.5, 0.015, 0.015)

# Run TISCA
tisca_results <- run_tisca(
  sim_func = run_simulation,
  comparison_pairs = comparison_pairs,
  mdes = mdes,
  target_power = 0.80,
  alpha = 0.05,
  batch_size = 100,
  initial_count = 100,
  correction_method = "holm",
  verbose = TRUE,
  save_results = TRUE
)

# Print final message
cat("\nTISCA analysis completed.\n")
cat("Required number of simulations:", tisca_results$J_final, "\n")
cat("Results saved to CSV files and plots generated.\n")
