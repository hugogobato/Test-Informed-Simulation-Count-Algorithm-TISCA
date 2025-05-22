#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TISCA: Test-Informed Simulation Count Algorithm
Generalized Python Package for Determining Optimal Simulation Counts

This is a Python implementation of the TISCA algorithm, translated from the original R code. This code is still in Beta version and any feedback is welcome
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from typing import Callable, List, Tuple, Dict, Union, Optional, Any


def welch_t_test(data: pd.DataFrame, metric_p: str, metric_b: str) -> Dict[str, float]:
    """
    Perform Welch's t-test for two independent samples.
    
    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing the metrics for both models
    metric_p : str
        Column name for the proposed model's metric
    metric_b : str
        Column name for the benchmark model's metric
    
    Returns
    -------
    dict
        Dictionary containing test results:
        - p_value: p-value from the test
        - t_stat: t-statistic
        - sd_p: standard deviation of proposed model's metric
        - sd_b: standard deviation of benchmark model's metric
        - mean_p: mean of proposed model's metric
        - mean_b: mean of benchmark model's metric
    """
    # Extract the metrics for both models
    values_p = data[metric_p].values
    values_b = data[metric_b].values
    
    # Perform Welch's t-test
    t_stat, p_value = stats.ttest_ind(values_p, values_b, equal_var=False)
    
    return {
        'p_value': p_value,
        't_stat': t_stat,
        'sd_p': np.std(values_p, ddof=1),
        'sd_b': np.std(values_b, ddof=1),
        'mean_p': np.mean(values_p),
        'mean_b': np.mean(values_b)
    }


def adjust_p_values(p_values: Union[List[float], np.ndarray], method: str) -> np.ndarray:
    """
    Adjust p-values for multiple testing using statsmodels.

    Parameters
    ----------
    p_values : list or numpy.ndarray
        List or array of p-values to adjust.
    method : str
        Method for adjustment. Supported methods by this wrapper:
        "none", "bonferroni", "holm", "fdr_bh" (Benjamini-Hochberg).
        For other methods supported by statsmodels, you can extend this function
        or call statsmodels.stats.multitest.multipletests directly.

    Returns
    -------
    numpy.ndarray
        Adjusted p-values.

    Raises
    ------
    ValueError
        If an unsupported correction method is provided.
    ImportError
        If statsmodels is not installed.
    """
    # Ensure p_values is a numpy array for statsmodels
    p_values_arr = np.asarray(p_values)

    if method == "none":
        return p_values_arr
    else:
        # These are the method strings directly accepted by statsmodels.stats.multitest.multipletests
        # that we are explicitly supporting in this wrapper.
        supported_methods_statsmodels = {
            "bonferroni",
            "holm",
            "fdr_bh"  # Benjamini-Hochberg
            # Add other statsmodels methods here if you want to expose them, e.g.:
            # "sidak", "holm-sidak", "fdr_by" (Benjamini-Yekutieli)
        }

        if method not in supported_methods_statsmodels:
            raise ValueError(
                f"Invalid or unsupported correction method: '{method}'. "
                f"This wrapper supports: 'none', {', '.join(sorted(list(supported_methods_statsmodels)))}"
            )

        # The `alpha` parameter in `multipletests` is used for determining
        # the `reject` boolean array. It does not affect the calculation
        # of `pvals_corrected` for these common methods.
        # We are primarily interested in the adjusted p-values here.
        reject, p_adjusted, _, _ = multipletests(
            p_values_arr,
            alpha=0.05,  # A common default, mainly for 'reject' output
            method=method # Pass the method string directly
        )

        return p_adjusted
    
def estimate_welch_power(n1: int, n2: int, sd1: float, sd2: float, delta: float, alpha: float) -> float:
    """
    Estimate statistical power for Welch's t-test.
    
    Parameters
    ----------
    n1 : int
        Sample size for first group
    n2 : int
        Sample size for second group
    sd1 : float
        Standard deviation for first group
    sd2 : float
        Standard deviation for second group
    delta : float
        Effect size (difference in means)
    alpha : float
        Significance level
    
    Returns
    -------
    float
        Estimated power
    """
    # Calculate degrees of freedom using Welch-Satterthwaite equation
    df = ((sd1**2/n1 + sd2**2/n2)**2) / \
         ((sd1**2/n1)**2/(n1-1) + (sd2**2/n2)**2/(n2-1))
    
    # Calculate non-centrality parameter
    ncp = delta / np.sqrt(sd1**2/n1 + sd2**2/n2)
    
    # Calculate critical value
    crit = stats.t.ppf(1-alpha/2, df)
    
    # Calculate power
    power = 1 - stats.t.cdf(crit, df, ncp) + stats.t.cdf(-crit, df, ncp)
    
    return power


def create_power_plot(power_tracking: pd.DataFrame, target_power: float, output_file: Optional[str] = None) -> plt.Figure:
    """
    Create a plot of power vs number of simulations.
    
    Parameters
    ----------
    power_tracking : pandas.DataFrame
        DataFrame with columns: simulations, comparison, power
    target_power : float
        Target power level to highlight
    output_file : str, optional
        File path to save the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    plt.figure(figsize=(10, 6))
    
    # Set the style
    sns.set_style("whitegrid")
    
    # Create the line plot
    ax = sns.lineplot(
        data=power_tracking, 
        x="simulations", 
        y="power", 
        hue="comparison",
        marker="o",
        linewidth=2
    )
    
    # Add target power line
    plt.axhline(y=target_power, linestyle="--", color="red", label=f"Target Power = {target_power}")
    
    # Set labels and title
    plt.title("Power vs Number of Simulations", fontsize=16, fontweight="bold")
    plt.xlabel("Number of Simulations", fontsize=12)
    plt.ylabel("Estimated Power", fontsize=12)
    plt.legend(title="Comparison", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Improve layout
    plt.tight_layout()
    
    # Save the plot if output_file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    
    return plt.gcf()


def create_pvalue_plot(pvalue_tracking: pd.DataFrame, alpha: float, J_final: int, 
                       output_file: Optional[str] = None) -> plt.Figure:
    """
    Create a bar plot of p-values for each comparison.
    
    Parameters
    ----------
    pvalue_tracking : pandas.DataFrame
        DataFrame with columns: simulations, comparison, p_raw, p_adj
    alpha : float
        Significance level to highlight
    J_final : int
        Final number of simulations
    output_file : str, optional
        File path to save the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    # Get the final p-values
    final_pvalues = pvalue_tracking[pvalue_tracking['simulations'] == pvalue_tracking['simulations'].max()].copy()
    
    # Transform p-values to -log10 scale for better visualization
    final_pvalues['neg_log10_p_raw'] = -np.log10(final_pvalues['p_raw'])
    final_pvalues['neg_log10_p_adj'] = -np.log10(final_pvalues['p_adj'])
    
    # Reshape for plotting
    plot_data = []
    for _, row in final_pvalues.iterrows():
        plot_data.append({
            'comparison': row['comparison'],
            'p_type': 'Raw p-value',
            'neg_log10_p': row['neg_log10_p_raw']
        })
        plot_data.append({
            'comparison': row['comparison'],
            'p_type': 'Adjusted p-value',
            'neg_log10_p': row['neg_log10_p_adj']
        })
    
    pvalue_plot_data = pd.DataFrame(plot_data)
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Set the style
    sns.set_style("whitegrid")
    
    # Create the bar plot
    ax = sns.barplot(
        data=pvalue_plot_data,
        x="comparison",
        y="neg_log10_p",
        hue="p_type",
        palette="Set2"
    )
    
    # Add significance line
    plt.axhline(y=-np.log10(alpha), linestyle="--", color="red")
    plt.text(
        0, 
        -np.log10(alpha) + 0.2, 
        f"α = {alpha}", 
        color="red", 
        ha="left", 
        va="bottom"
    )
    
    # Set labels and title
    plt.title("Significance of Comparisons", fontsize=16, fontweight="bold")
    plt.suptitle(f"After {J_final} Simulations", fontsize=12)
    plt.xlabel("Comparison", fontsize=12)
    plt.ylabel("-log10(p-value)", fontsize=12)
    plt.legend(title="P-value Type")
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    
    # Improve layout
    plt.tight_layout()
    
    # Save the plot if output_file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    
    return plt.gcf()


def run_tisca(
    sim_func: Callable[[int], pd.DataFrame],  # Function that runs one simulation and returns metrics
    comparison_pairs: List[Tuple[str, str]],   # List of pairs to compare (metric_p, metric_b)
    mdes: List[float],                         # List of minimum detectable effect sizes
    target_power: float = 0.80,                # Target statistical power
    alpha: float = 0.05,                       # Significance level
    batch_size: int = 50,                      # Number of simulations per batch
    initial_count: Optional[int] = None,       # Initial simulation count (defaults to batch_size if None)
    correction_method: str = "none",           # Multiple testing correction method
    verbose: bool = True,                      # Whether to print progress updates
    save_results: bool = True,                 # Whether to save results to CSV files
    output_dir: str = "."                      # Directory to save output files
) -> Dict[str, Any]:
    """
    Python implementation of the Test-Informed Simulation Count Algorithm (TISCA)
    for determining the optimal number of simulations required to achieve a target
    statistical power when comparing model performance metrics.
    
    Parameters
    ----------
    sim_func : callable
        A function that takes a seed parameter and returns a pandas DataFrame
        containing performance metrics for all models being compared in a single simulation.
        The function must have the signature: function(seed) -> DataFrame
    comparison_pairs : list of tuples
        A list of pairs, where each pair is a tuple of two strings:
        ("metric_proposed", "metric_benchmark"). These strings must match column names
        in the DataFrame returned by sim_func.
    mdes : list of float
        A list of Minimum Detectable Effect Sizes, one for each comparison pair.
        For metrics where lower is better (e.g., RMSE, PEHE), use negative values.
        For metrics where higher is better (e.g., coverage), use positive values.
    target_power : float, optional
        The desired statistical power (1-β), typically 0.80. Default is 0.80.
    alpha : float, optional
        The significance level (α), typically 0.05. Default is 0.05.
    batch_size : int, optional
        The number of simulations to run in each iteration. Default is 50.
    initial_count : int, optional
        The number of simulations to run before the first power check.
        Defaults to batch_size if None.
    correction_method : str, optional
        The method for adjusting p-values in multiple testing.
        Options: "none", "bonferroni", "holm", "BH" (Benjamini-Hochberg).
        Default is "none".
    verbose : bool, optional
        Whether to print progress updates. Default is True.
    save_results : bool, optional
        Whether to save results to CSV files. Default is True.
    output_dir : str, optional
        Directory to save output files. Default is current directory.
    
    Returns
    -------
    dict
        A dictionary containing:
        - J_final: The final number of simulations required
        - P_raw: List of raw p-values
        - P_adj: List of adjusted p-values
        - Stats: List of test statistics
        - P_achieved: List of achieved powers
        - results: DataFrame with all simulation results
        - comparison_pairs: The comparison pairs used
        - mdes: The minimum detectable effect sizes used
        - summary_table: A DataFrame summarizing the results
        - power_tracking: DataFrame tracking power estimates across simulation counts
        - pvalue_tracking: DataFrame tracking p-values across simulation counts
    """
    # Validate inputs
    if not callable(sim_func):
        raise TypeError("sim_func must be a function")
    
    if not isinstance(comparison_pairs, list) or len(comparison_pairs) == 0:
        raise ValueError("comparison_pairs must be a non-empty list")
    
    for pair in comparison_pairs:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError("Each comparison pair must be a tuple of length 2")
    
    if not isinstance(mdes, list) or len(mdes) != len(comparison_pairs):
        raise ValueError("mdes must be a list with the same length as comparison_pairs")
    
    if not isinstance(target_power, (int, float)) or target_power <= 0 or target_power >= 1:
        raise ValueError("target_power must be a number between 0 and 1")
    
    if not isinstance(alpha, (int, float)) or alpha <= 0 or alpha >= 1:
        raise ValueError("alpha must be a number between 0 and 1")
    
    if not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    
    if initial_count is not None and (not isinstance(initial_count, int) or initial_count <= 0):
        raise ValueError("initial_count must be a positive integer or None")
    
    if correction_method not in ["none", "bonferroni", "holm", "BH"]:
        raise ValueError("correction_method must be one of: 'none', 'bonferroni', 'holm', 'BH'")
    
    # Set initial count to batch_size if not specified
    if initial_count is None:
        initial_count = batch_size
    
    # Create output directory if it doesn't exist and save_results is True
    if save_results and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Initialization Phase
    J = 0
    results_agg = pd.DataFrame()  # Empty dataframe to store results
    K = len(comparison_pairs)
    P_current = np.zeros(K)       # Array of K zeros for current power
    
    # Create data structures to track power and p-values over iterations
    power_tracking = pd.DataFrame(columns=["simulations", "comparison", "power"])
    pvalue_tracking = pd.DataFrame(columns=["simulations", "comparison", "p_raw", "p_adj"])
    
    # Print initial information if verbose
    if verbose:
        print("Starting TISCA algorithm...")
        print(f"Target power: {target_power}")
        print(f"Significance level: {alpha}")
        print(f"Batch size: {batch_size}")
        print(f"Initial count: {initial_count}")
        print(f"Correction method: {correction_method}")
        print(f"Number of comparisons: {K}\n")
    
    # Optional Initial Run
    if initial_count > 0 and initial_count >= batch_size:
        if verbose:
            print(f"Running initial {initial_count} simulations...")
            progress_bar = tqdm(total=initial_count, disable=not verbose)
        
        # Run initial simulations
        for i in range(1, initial_count + 1):
            if verbose:
                progress_bar.update(1)
            
            # Try to run the simulation with error handling
            try:
                metrics_run = sim_func(seed=i)
                
                # Validate that the simulation function returns a DataFrame
                if not isinstance(metrics_run, pd.DataFrame):
                    raise TypeError("Simulation function must return a pandas DataFrame")
                
                # Validate that all required metrics are in the DataFrame
                all_metrics = set([metric for pair in comparison_pairs for metric in pair])
                missing_metrics = [metric for metric in all_metrics if metric not in metrics_run.columns]
                if missing_metrics:
                    raise ValueError(f"Missing metrics in simulation output: {', '.join(missing_metrics)}")
                
                # Add to results
                results_agg = pd.concat([results_agg, metrics_run], ignore_index=True)
            
            except Exception as e:
                if verbose:
                    progress_bar.close()
                raise RuntimeError(f"Error in simulation function: {str(e)}")
        
        if verbose:
            progress_bar.close()
        
        J = initial_count
        
        # Perform initial power calculation
        p_values_raw = np.zeros(K)
        test_stats = np.zeros(K)
        mean_diffs = np.zeros(K)
        
        for k in range(K):
            # Extract comparison information
            metric_p, metric_b = comparison_pairs[k]
            comparison_name = f"{metric_p} vs {metric_b}"
            
            # Perform Welch's t-test
            test_result = welch_t_test(results_agg, metric_p, metric_b)
            
            # Store results
            p_values_raw[k] = test_result['p_value']
            test_stats[k] = test_result['t_stat']
            mean_diffs[k] = test_result['mean_p'] - test_result['mean_b']
            
            # Estimate power if standard deviations are valid
            if test_result['sd_p'] > 0 and test_result['sd_b'] > 0:
                P_current[k] = estimate_welch_power(
                    J, J, 
                    test_result['sd_p'], test_result['sd_b'], 
                    mdes[k], alpha
                )
            else:
                P_current[k] = 0
            
            # Track power
            power_tracking = pd.concat([
                power_tracking,
                pd.DataFrame({
                    "simulations": [J],
                    "comparison": [comparison_name],
                    "power": [P_current[k]]
                })
            ], ignore_index=True)
        
        # Adjust p-values for multiple testing
        p_values_adj = adjust_p_values(p_values_raw, correction_method)
        
        # Track p-values
        for k in range(K):
            comparison_name = f"{comparison_pairs[k][0]} vs {comparison_pairs[k][1]}"
            pvalue_tracking = pd.concat([
                pvalue_tracking,
                pd.DataFrame({
                    "simulations": [J],
                    "comparison": [comparison_name],
                    "p_raw": [p_values_raw[k]],
                    "p_adj": [p_values_adj[k]]
                })
            ], ignore_index=True)
        
        if verbose:
            print(f"\nAfter initial {J} simulations:")
            for k in range(K):
                print(f"Comparison {k+1}: Power = {P_current[k]:.4f}, "
                      f"p-value = {p_values_raw[k]:.4f}, "
                      f"adjusted p-value = {p_values_adj[k]:.4f}")
            print(f"Minimum power: {np.min(P_current):.4f}\n")
    
    # Iterative Simulation and Power Check Loop
    iteration = 1
    
    while np.min(P_current) < target_power:
        if verbose:
            print(f"Iteration {iteration}: Running batch of {batch_size} simulations "
                  f"(total will be {J + batch_size})...")
            progress_bar = tqdm(total=batch_size, disable=not verbose)
        
        # Run a new batch of simulations
        start_seed = J + 1
        end_seed = J + batch_size
        
        for i in range(start_seed, end_seed + 1):
            if verbose:
                progress_bar.update(1)
            
            # Try to run the simulation with error handling
            try:
                metrics_run = sim_func(seed=i)
                
                # Validate that the simulation function returns a DataFrame
                if not isinstance(metrics_run, pd.DataFrame):
                    raise TypeError("Simulation function must return a pandas DataFrame")
                
                # Validate that all required metrics are in the DataFrame
                all_metrics = set([metric for pair in comparison_pairs for metric in pair])
                missing_metrics = [metric for metric in all_metrics if metric not in metrics_run.columns]
                if missing_metrics:
                    raise ValueError(f"Missing metrics in simulation output: {', '.join(missing_metrics)}")
                
                # Add to results
                results_agg = pd.concat([results_agg, metrics_run], ignore_index=True)
            
            except Exception as e:
                if verbose:
                    progress_bar.close()
                raise RuntimeError(f"Error in simulation function: {str(e)}")
        
        if verbose:
            progress_bar.close()
        
        J += batch_size
        
        # Perform tests and estimate power on current J simulations
        p_values_raw = np.zeros(K)
        test_stats = np.zeros(K)
        mean_diffs = np.zeros(K)
        
        for k in range(K):
            # Extract comparison information
            metric_p, metric_b = comparison_pairs[k]
            comparison_name = f"{metric_p} vs {metric_b}"
            
            # Perform Welch's t-test
            test_result = welch_t_test(results_agg, metric_p, metric_b)
            
            # Store results
            p_values_raw[k] = test_result['p_value']
            test_stats[k] = test_result['t_stat']
            mean_diffs[k] = test_result['mean_p'] - test_result['mean_b']
            
            # Estimate power if standard deviations are valid
            if test_result['sd_p'] > 0 and test_result['sd_b'] > 0:
                P_current[k] = estimate_welch_power(
                    J, J, 
                    test_result['sd_p'], test_result['sd_b'], 
                    mdes[k], alpha
                )
            else:
                P_current[k] = 0
            
            # Track power
            power_tracking = pd.concat([
                power_tracking,
                pd.DataFrame({
                    "simulations": [J],
                    "comparison": [comparison_name],
                    "power": [P_current[k]]
                })
            ], ignore_index=True)
        
        # Adjust p-values for multiple testing
        p_values_adj = adjust_p_values(p_values_raw, correction_method)
        
        # Track p-values
        for k in range(K):
            comparison_name = f"{comparison_pairs[k][0]} vs {comparison_pairs[k][1]}"
            pvalue_tracking = pd.concat([
                pvalue_tracking,
                pd.DataFrame({
                    "simulations": [J],
                    "comparison": [comparison_name],
                    "p_raw": [p_values_raw[k]],
                    "p_adj": [p_values_adj[k]]
                })
            ], ignore_index=True)
        
        if verbose:
            print(f"\nAfter {J} simulations:")
            for k in range(K):
                print(f"Comparison {k+1}: Power = {P_current[k]:.4f}, "
                      f"p-value = {p_values_raw[k]:.4f}, "
                      f"adjusted p-value = {p_values_adj[k]:.4f}")
            print(f"Minimum power: {np.min(P_current):.4f}\n")
        
        iteration += 1
    
    # Output Phase
    J_final = J
    P_raw = p_values_raw
    P_adj = p_values_adj
    Stats = test_stats
    P_achieved = P_current
    
    # Create summary table
    summary_table = pd.DataFrame(columns=[
        "Comparison", "Proposed_Metric", "Benchmark_Metric", "MDE", "Mean_Difference",
        "Power", "P_Value", "Adjusted_P_Value", "T_Statistic"
    ])
    
    for k in range(len(comparison_pairs)):
        summary_table.loc[k] = [
            f"{comparison_pairs[k][0]} vs {comparison_pairs[k][1]}",
            comparison_pairs[k][0],
            comparison_pairs[k][1],
            mdes[k],
            mean_diffs[k],
            P_achieved[k],
            P_raw[k],
            P_adj[k],
            Stats[k]
        ]
    
    if verbose:
        print("\n=== TISCA COMPLETED ===")
        print(f"Final number of simulations required: {J_final}")
        print("Final achieved powers:")
        for k in range(K):
            print(f"Comparison {k+1} ({comparison_pairs[k][0]} vs {comparison_pairs[k][1]}): "
                  f"Power = {P_achieved[k]:.4f}, p-value = {P_raw[k]:.4f}, "
                  f"adjusted p-value = {P_adj[k]:.4f}")
    
    # Create visualizations and save results if requested
    if save_results:
        # Ensure output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # File paths
        results_file = os.path.join(output_dir, "tisca_simulation_results.csv")
        summary_file = os.path.join(output_dir, "tisca_summary_results.csv")
        power_tracking_file = os.path.join(output_dir, "power_tracking.csv")
        pvalue_tracking_file = os.path.join(output_dir, "pvalue_tracking.csv")
        power_plot_file = os.path.join(output_dir, "power_vs_simulations.png")
        pvalue_plot_file = os.path.join(output_dir, "pvalue_comparison.png")
        
        # Save data to CSV files
        results_agg.to_csv(results_file, index=False)
        summary_table.to_csv(summary_file, index=False)
        power_tracking.to_csv(power_tracking_file, index=False)
        pvalue_tracking.to_csv(pvalue_tracking_file, index=False)
        
        # Create and save visualizations
        try:
            # Plot power vs simulations
            create_power_plot(power_tracking, target_power, power_plot_file)
            
            # Plot p-values
            create_pvalue_plot(pvalue_tracking, alpha, J_final, pvalue_plot_file)
            
            if verbose:
                print("\nResults saved to CSV files:")
                print(f"- {results_file}")
                print(f"- {summary_file}")
                print(f"- {power_tracking_file}")
                print(f"- {pvalue_tracking_file}")
                print("\nVisualizations saved:")
                print(f"- {power_plot_file}")
                print(f"- {pvalue_plot_file}")
        except Exception as e:
            if verbose:
                print(f"\nResults saved to CSV files, but visualizations could not be created: {str(e)}")
    
    # Return results
    return {
        "J_final": J_final,
        "P_raw": P_raw,
        "P_adj": P_adj,
        "Stats": Stats,
        "P_achieved": P_achieved,
        "results": results_agg,
        "comparison_pairs": comparison_pairs,
        "mdes": mdes,
        "summary_table": summary_table,
        "power_tracking": power_tracking,
        "pvalue_tracking": pvalue_tracking
    }


# Example usage
if __name__ == "__main__":
    # This code will only run if the script is executed directly, not when imported
    
    # Define a simple simulation function
    def sim_func(seed):
        np.random.seed(seed)
        # Simulate data and fit models
        # Return performance metrics as a DataFrame
        return pd.DataFrame({
            'model_a_rmse': [np.random.normal(0.5, 0.1)],
            'model_b_rmse': [np.random.normal(0.6, 0.1)],
            'model_a_coverage': [np.random.normal(0.9, 0.05)],
            'model_b_coverage': [np.random.normal(0.85, 0.05)]
        })
    
    # Define comparison pairs
    comparison_pairs = [
        ('model_a_rmse', 'model_b_rmse'),
        ('model_a_coverage', 'model_b_coverage')
    ]
    
    # Define minimum detectable effect sizes
    # Negative for RMSE (lower is better), positive for coverage (higher is better)
    mdes = [-0.1, 0.05]
    
    # Run TISCA
    results = run_tisca(
        sim_func=sim_func,
        comparison_pairs=comparison_pairs,
        mdes=mdes,
        target_power=0.80,
        alpha=0.05,
        batch_size=10,
        initial_count=20,
        verbose=True,
        save_results=True
    )
    
    # Access results
    print(f"\nRequired number of simulations: {results['J_final']}")
    print(results['summary_table'])
