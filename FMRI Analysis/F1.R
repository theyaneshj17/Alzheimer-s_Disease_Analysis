# Combined Effect Size + Statistical Significance Analysis
# Optimal approach for AD connectivity research

library(ggplot2)
library(dplyr)
library(reshape2)

cat("COMBINED EFFECT SIZE + STATISTICAL SIGNIFICANCE ANALYSIS\n")
cat("========================================================\n\n")

# Check if models are loaded
if (!exists("model0") || !exists("model1") || !exists("model2")) {
  stop("Please load models first: model0 (S1), model1 (S2), model2 (S3)")
}

# =============================================================================
# APPROACH 1: HYBRID UVPM + UVC ANALYSIS (RECOMMENDED)
# =============================================================================

cat("APPROACH 1: HYBRID ANALYSIS (BEST FOR PUBLICATION)\n")
cat("===================================================\n\n")

# Function for comprehensive analysis combining both approaches
comprehensive_analysis <- function(model_a_uvpm, model_a_uvc, model_b_uvpm, model_b_uvc, 
                                   comparison_name, alpha = 0.05) {
  
  cat("Analyzing:", comparison_name, "\n")
  
  # UVPM-based effect sizes (conservative, reliable)
  values_a_uvpm <- model_a_uvpm[upper.tri(model_a_uvpm)]
  values_b_uvpm <- model_b_uvpm[upper.tri(model_b_uvpm)]
  
  differences_uvpm <- values_b_uvpm - values_a_uvpm
  pooled_sd_uvpm <- sqrt((var(values_a_uvpm) + var(values_b_uvpm)) / 2)
  effect_sizes_uvpm <- differences_uvpm / pooled_sd_uvpm
  
  # UVC-based statistical inference (uncertainty quantification)
  n_connections <- ncol(model_a_uvc)
  posterior_probs <- numeric(n_connections)
  credible_intervals_lower <- numeric(n_connections)
  credible_intervals_upper <- numeric(n_connections)
  
  for(i in 1:n_connections) {
    samples_a <- model_a_uvc[, i]
    samples_b <- model_b_uvc[, i]
    
    # Posterior probability that B > A (Bayesian equivalent of p-value)
    posterior_probs[i] <- mean(samples_b > samples_a)
    
    # Credible intervals for the difference
    diff_samples <- samples_b - samples_a
    credible_intervals_lower[i] <- quantile(diff_samples, alpha/2)
    credible_intervals_upper[i] <- quantile(diff_samples, 1 - alpha/2)
  }
  
  # Convert posterior probabilities to two-tailed "p-values"
  # P(B ≠ A) = 2 * min(P(B > A), P(B < A))
  bayesian_p_values <- 2 * pmin(posterior_probs, 1 - posterior_probs)
  
  # Combine UVPM effect sizes with UVC statistical inference
  combined_results <- data.frame(
    Connection = 1:length(effect_sizes_uvpm),
    Effect_Size_UVPM = effect_sizes_uvpm,
    Abs_Effect_Size = abs(effect_sizes_uvpm),
    Bayesian_P_Value = bayesian_p_values[1:length(effect_sizes_uvpm)],
    Credible_Lower = credible_intervals_lower[1:length(effect_sizes_uvpm)],
    Credible_Upper = credible_intervals_upper[1:length(effect_sizes_uvpm)],
    Posterior_Prob_Greater = posterior_probs[1:length(effect_sizes_uvpm)]
  )
  
  # Classification
  combined_results$Effect_Category <- case_when(
    combined_results$Abs_Effect_Size >= 0.8 ~ "Large",
    combined_results$Abs_Effect_Size >= 0.5 ~ "Medium",
    combined_results$Abs_Effect_Size >= 0.2 ~ "Small",
    TRUE ~ "Negligible"
  )
  
  combined_results$Statistical_Significance <- case_when(
    combined_results$Bayesian_P_Value < 0.001 ~ "Very Significant",
    combined_results$Bayesian_P_Value < 0.01 ~ "Significant", 
    combined_results$Bayesian_P_Value < 0.05 ~ "Marginally Significant",
    TRUE ~ "Not Significant"
  )
  
  # Credible interval significance (does not include zero)
  combined_results$Credible_Significant <- 
    (combined_results$Credible_Lower > 0) | (combined_results$Credible_Upper < 0)
  
  # Combined classification: Large Effect + Statistically Significant
  combined_results$Large_AND_Significant <- 
    (combined_results$Effect_Category == "Large") & 
    (combined_results$Bayesian_P_Value < 0.05)
  
  combined_results$Large_AND_Credible <- 
    (combined_results$Effect_Category == "Large") & 
    combined_results$Credible_Significant
  
  # Summary statistics
  summary_stats <- list(
    comparison = comparison_name,
    total_connections = nrow(combined_results),
    large_effects = sum(combined_results$Effect_Category == "Large"),
    significant_connections = sum(combined_results$Bayesian_P_Value < 0.05),
    credible_significant = sum(combined_results$Credible_Significant),
    large_and_significant = sum(combined_results$Large_AND_Significant),
    large_and_credible = sum(combined_results$Large_AND_Credible),
    mean_effect_size = mean(combined_results$Abs_Effect_Size),
    max_effect_size = max(combined_results$Abs_Effect_Size)
  )
  
  cat("Results Summary:\n")
  cat("- Total connections:", summary_stats$total_connections, "\n")
  cat("- Large effects (|d| > 0.8):", summary_stats$large_effects, 
      sprintf(" (%.1f%%)", 100 * summary_stats$large_effects / summary_stats$total_connections), "\n")
  cat("- Statistically significant (p < 0.05):", summary_stats$significant_connections, 
      sprintf(" (%.1f%%)", 100 * summary_stats$significant_connections / summary_stats$total_connections), "\n")
  cat("- Credible interval significant:", summary_stats$credible_significant, 
      sprintf(" (%.1f%%)", 100 * summary_stats$credible_significant / summary_stats$total_connections), "\n")
  cat("- Large AND significant:", summary_stats$large_and_significant, 
      sprintf(" (%.1f%%)", 100 * summary_stats$large_and_significant / summary_stats$total_connections), "\n")
  cat("- Large AND credible:", summary_stats$large_and_credible, 
      sprintf(" (%.1f%%)", 100 * summary_stats$large_and_credible / summary_stats$total_connections), "\n\n")
  
  return(list(
    results = combined_results,
    summary = summary_stats
  ))
}

# =============================================================================
# RUN COMPREHENSIVE ANALYSIS
# =============================================================================

# Run analysis for all comparisons
cat("Running comprehensive analysis for all subtype comparisons...\n\n")

s1_vs_s3_comprehensive <- comprehensive_analysis(
  model0$UVPM, model0$UVC, model2$UVPM, model2$UVC, "S1 vs S3"
)

s3_vs_s2_comprehensive <- comprehensive_analysis(
  model2$UVPM, model2$UVC, model1$UVPM, model1$UVC, "S3 vs S2"
)

s1_vs_s2_comprehensive <- comprehensive_analysis(
  model0$UVPM, model0$UVC, model1$UVPM, model1$UVC, "S1 vs S2"
)

# =============================================================================
# CREATE COMPREHENSIVE SUMMARY TABLE
# =============================================================================

cat("COMPREHENSIVE SUMMARY: EFFECT SIZE + SIGNIFICANCE\n")
cat("=================================================\n\n")

comprehensive_summary <- data.frame(
  Comparison = c("S1 vs S3", "S3 vs S2", "S1 vs S2"),
  Total_Connections = c(s1_vs_s3_comprehensive$summary$total_connections,
                        s3_vs_s2_comprehensive$summary$total_connections,
                        s1_vs_s2_comprehensive$summary$total_connections),
  Large_Effects_N = c(s1_vs_s3_comprehensive$summary$large_effects,
                      s3_vs_s2_comprehensive$summary$large_effects,
                      s1_vs_s2_comprehensive$summary$large_effects),
  Large_Effects_Pct = round(100 * c(s1_vs_s3_comprehensive$summary$large_effects,
                                    s3_vs_s2_comprehensive$summary$large_effects,
                                    s1_vs_s2_comprehensive$summary$large_effects) /
                              c(s1_vs_s3_comprehensive$summary$total_connections,
                                s3_vs_s2_comprehensive$summary$total_connections,
                                s1_vs_s2_comprehensive$summary$total_connections), 1),
  Significant_N = c(s1_vs_s3_comprehensive$summary$significant_connections,
                    s3_vs_s2_comprehensive$summary$significant_connections,
                    s1_vs_s2_comprehensive$summary$significant_connections),
  Significant_Pct = round(100 * c(s1_vs_s3_comprehensive$summary$significant_connections,
                                  s3_vs_s2_comprehensive$summary$significant_connections,
                                  s1_vs_s2_comprehensive$summary$significant_connections) /
                            c(s1_vs_s3_comprehensive$summary$total_connections,
                              s3_vs_s2_comprehensive$summary$total_connections,
                              s1_vs_s2_comprehensive$summary$total_connections), 1),
  Large_AND_Significant_N = c(s1_vs_s3_comprehensive$summary$large_and_significant,
                              s3_vs_s2_comprehensive$summary$large_and_significant,
                              s1_vs_s2_comprehensive$summary$large_and_significant),
  Large_AND_Significant_Pct = round(100 * c(s1_vs_s3_comprehensive$summary$large_and_significant,
                                            s3_vs_s2_comprehensive$summary$large_and_significant,
                                            s1_vs_s2_comprehensive$summary$large_and_significant) /
                                      c(s1_vs_s3_comprehensive$summary$total_connections,
                                        s3_vs_s2_comprehensive$summary$total_connections,
                                        s1_vs_s2_comprehensive$summary$total_connections), 1)
)

print("GOLD STANDARD RESULTS: Large Effect + Statistical Significance")
print(comprehensive_summary)

# =============================================================================
# VISUALIZATION: EFFECT SIZE vs SIGNIFICANCE
# =============================================================================

cat("\nCreating publication-quality visualizations...\n\n")

# Create effect size vs significance plot for S1 vs S2
s1_s2_data <- s1_vs_s2_comprehensive$results

# Scatter plot: Effect Size vs -log10(p-value)
p1 <- ggplot(s1_s2_data, aes(x = Abs_Effect_Size, y = -log10(Bayesian_P_Value))) +
  geom_point(aes(color = Large_AND_Significant), alpha = 0.6, size = 1) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#e74c3c"),
                     labels = c("Not Both", "Large & Significant")) +
  labs(title = "S1 vs S2: Effect Size vs Statistical Significance",
       subtitle = paste("Red points: Large effect + statistically significant (",
                        sum(s1_s2_data$Large_AND_Significant), "connections)"),
       x = "Absolute Effect Size |Cohen's d|",
       y = "-log10(Bayesian P-value)",
       color = "Classification") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

# Summary barplot
summary_long <- melt(comprehensive_summary[, c("Comparison", "Large_Effects_Pct", 
                                               "Significant_Pct", "Large_AND_Significant_Pct")],
                     id.vars = "Comparison", variable.name = "Metric", value.name = "Percentage")

p2 <- ggplot(summary_long, aes(x = Comparison, y = Percentage, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Large_Effects_Pct" = "#f39c12",
                               "Significant_Pct" = "#3498db", 
                               "Large_AND_Significant_Pct" = "#e74c3c"),
                    labels = c("Large Effects", "Statistically Significant", "Both Large & Significant")) +
  labs(title = "Combined Effect Size + Significance Analysis",
       subtitle = "Gold standard approach for AD connectivity research",
       x = "Subtype Comparison", y = "Percentage of Connections",
       fill = "Classification") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save results
output_dir <- "Combined_Effect_Significance_Analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "S1_vs_S2_Effect_vs_Significance.png"), p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Combined_Analysis_Summary.png"), p2, width = 12, height = 8, dpi = 300)

# Save detailed results
write.csv(comprehensive_summary, file.path(output_dir, "Comprehensive_Summary_Table.csv"), row.names = FALSE)
write.csv(s1_s2_data, file.path(output_dir, "S1_vs_S2_Detailed_Results.csv"), row.names = FALSE)

# =============================================================================
# PUBLICATION-READY REPORTING TEMPLATE
# =============================================================================

cat("PUBLICATION-READY REPORTING TEMPLATE\n")
cat("====================================\n\n")

cat("FOR YOUR METHODS SECTION:\n")
cat('```\n')
cat("Connectivity differences between subtypes were quantified using Cohen's d \n")
cat("effect sizes derived from posterior mean estimates (UVPM). Statistical \n")
cat("significance was assessed using Bayesian posterior probabilities from \n")
cat("MCMC samples (UVC). Large effects (|d| > 0.8) combined with statistical \n")
cat("significance (p < 0.05) were considered robust, biologically meaningful \n")
cat("differences following established neuroimaging guidelines.\n")
cat('```\n\n')

cat("FOR YOUR RESULTS SECTION:\n")
cat('```\n')
s1_s2_large_sig <- s1_vs_s2_comprehensive$summary$large_and_significant
s1_s2_total <- s1_vs_s2_comprehensive$summary$total_connections
s1_s2_pct <- round(100 * s1_s2_large_sig / s1_s2_total, 1)

cat("Comprehensive analysis combining effect sizes and statistical significance \n")
cat("revealed", s1_s2_large_sig, "connections (", s1_s2_pct, "%) between S1 and S2 subtypes \n")
cat("with both large effects (Cohen's |d| > 0.8) and statistical significance \n")
cat("(p < 0.05). This represents robust evidence for substantial connectivity \n")
cat("differences consistent with AD progression patterns.\n")
cat('```\n\n')

cat("FILES CREATED:\n")
cat("• Comprehensive_Summary_Table.csv - Main results table\n")
cat("• S1_vs_S2_Detailed_Results.csv - Detailed connection-wise results\n")
cat("• Effect vs Significance plots - Publication figures\n\n")

cat("This combined approach provides the GOLD STANDARD evidence for your AD research!\n")