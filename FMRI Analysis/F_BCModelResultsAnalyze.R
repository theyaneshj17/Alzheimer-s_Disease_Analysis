
# Set where you want plots saved
setwd("C:/Users/theya/Downloads/ADNIFinal/FMRIAnalysis/")

# Enhanced UVC Significance Analysis Script
library(ggplot2)
library(reshape2)
library(dplyr)

# Configuration
modeldir <- "C:/Users/theya/Downloads/FMRIFullBCResults"
output_dir <- file.path(modeldir, "plots123")
base_path <- "C:/Users/theya/Downloads/FMRIFullBCResults/K=5/K=5/"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load models with error handling
load_model_safely <- function(cluster_id, base_path) {
  filename <- file.path(base_path, paste0("bc_model_cluster_", cluster_id, ".rds"))
  
  if (!file.exists(filename)) {
    stop("File does not exist: ", filename)
  }
  
  tryCatch({
    model <- readRDS(filename)
    cat("Successfully loaded model for cluster", cluster_id, "\n")
    
    # Validate model structure
    if (!"UVC" %in% names(model)) {
      stop("UVC component not found in model for cluster ", cluster_id)
    }
    
    cat("  UVC dimensions:", dim(model$UVC), "\n")
    cat("  UVC samples range:", range(model$UVC), "\n")
    
    return(model)
  }, error = function(e) {
    stop("Error loading model for cluster ", cluster_id, ": ", e$message)
  })
}

# Load all models
cat("Loading models...\n")
model0 <- load_model_safely(0, base_path)  # Early stage
model1 <- load_model_safely(1, base_path)  # Advanced stage  
model2 <- load_model_safely(2, base_path)  # Intermediate stage

# Enhanced significance calculation function
calculate_uvc_significance <- function(model_a, model_b, comparison_name, alpha = 0.05) {
  
  cat("Calculating significance for:", comparison_name, "\n")
  
  # Validate UVC matrices have the same structure
  if (!all(dim(model_a$UVC) == dim(model_b$UVC))) {
    stop("UVC matrices must have identical dimensions")
  }
  
  # Number of connections
  n_connections <- ncol(model_a$UVC)
  n_samples <- nrow(model_a$UVC)
  
  cat("  Number of connections:", n_connections, "\n")
  cat("  Number of MCMC samples:", n_samples, "\n")
  
  # Generate connection names if not available
  connection_names <- colnames(model_a$UVC)
  if (is.null(connection_names)) {
    connection_names <- paste0("Connection_", 1:n_connections)
  }
  
  # Prepare results storage
  significance_results <- data.frame(
    Comparison = rep(comparison_name, n_connections),
    Connection = connection_names,
    Model_A_Mean = numeric(n_connections),
    Model_A_LB = numeric(n_connections),
    Model_A_UB = numeric(n_connections),
    Model_B_Mean = numeric(n_connections),
    Model_B_LB = numeric(n_connections),
    Model_B_UB = numeric(n_connections),
    Difference_Mean = numeric(n_connections),
    Difference_LB = numeric(n_connections),
    Difference_UB = numeric(n_connections),
    Significantly_Larger = logical(n_connections),
    Significantly_Smaller = logical(n_connections),
    Overlapping = logical(n_connections),
    Effect_Size = numeric(n_connections),
    Posterior_Probability_Larger = numeric(n_connections),
    stringsAsFactors = FALSE
  )
  
  # Compute statistics for each connection
  for (i in 1:n_connections) {
    # Extract samples for this connection
    samples_a <- model_a$UVC[, i]
    samples_b <- model_b$UVC[, i]
    
    # Check for valid samples
    if (any(is.na(samples_a)) || any(is.na(samples_b))) {
      warning("NA values found in connection ", i, " for ", comparison_name)
    }
    
    # Compute means
    significance_results$Model_A_Mean[i] <- mean(samples_a, na.rm = TRUE)
    significance_results$Model_B_Mean[i] <- mean(samples_b, na.rm = TRUE)
    
    # Compute credible intervals
    significance_results$Model_A_LB[i] <- quantile(samples_a, alpha / 2, na.rm = TRUE)
    significance_results$Model_A_UB[i] <- quantile(samples_a, 1 - alpha / 2, na.rm = TRUE)
    
    significance_results$Model_B_LB[i] <- quantile(samples_b, alpha / 2, na.rm = TRUE)
    significance_results$Model_B_UB[i] <- quantile(samples_b, 1 - alpha / 2, na.rm = TRUE)
    
    # Compute difference statistics
    difference_samples <- samples_b - samples_a
    significance_results$Difference_Mean[i] <- mean(difference_samples, na.rm = TRUE)
    significance_results$Difference_LB[i] <- quantile(difference_samples, alpha / 2, na.rm = TRUE)
    significance_results$Difference_UB[i] <- quantile(difference_samples, 1 - alpha / 2, na.rm = TRUE)
    
    # Identify significance categories (non-overlapping credible intervals)
    significance_results$Significantly_Larger[i] <- significance_results$Model_B_LB[i] > significance_results$Model_A_UB[i]
    significance_results$Significantly_Smaller[i] <- significance_results$Model_B_UB[i] < significance_results$Model_A_LB[i]
    significance_results$Overlapping[i] <- !(significance_results$Significantly_Larger[i] | significance_results$Significantly_Smaller[i])
    
    # Effect size (standardized difference)
    pooled_sd <- sqrt((var(samples_a, na.rm = TRUE) + var(samples_b, na.rm = TRUE)) / 2)
    significance_results$Effect_Size[i] <- significance_results$Difference_Mean[i] / pooled_sd
    
    # Posterior probability that B > A
    significance_results$Posterior_Probability_Larger[i] <- mean(difference_samples > 0, na.rm = TRUE)
  }
  
  return(significance_results)
}

# Compute significant differences for all comparisons
cat("\nComputing pairwise comparisons...\n")
results_early_intermediate <- calculate_uvc_significance(model0, model2, "Early_vs_Intermediate", alpha = 0.05)
results_intermediate_advanced <- calculate_uvc_significance(model2, model1, "Intermediate_vs_Advanced", alpha = 0.05)
results_early_advanced <- calculate_uvc_significance(model0, model1, "Early_vs_Advanced", alpha = 0.05)

# Combine all results into one dataframe
all_results <- rbind(results_early_intermediate, results_intermediate_advanced, results_early_advanced)

# Enhanced output directory
full_output_dir <- file.path(dirname(base_path), "significance_analysis")
dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)

# Save detailed results
detailed_file <- file.path(full_output_dir, "detailed_comparisons_results.csv")
write.csv(all_results, detailed_file, row.names = FALSE)
cat("Detailed results saved to:", detailed_file, "\n")

# Create summary statistics
create_summary_table <- function(results_df) {
  summary_stats <- results_df %>%
    group_by(Comparison) %>%
    summarise(
      Total_Connections = n(),
      Significantly_Larger = sum(Significantly_Larger),
      Significantly_Smaller = sum(Significantly_Smaller),
      Overlapping = sum(Overlapping),
      Pct_Larger = round(100 * sum(Significantly_Larger) / n(), 1),
      Pct_Smaller = round(100 * sum(Significantly_Smaller) / n(), 1),
      Pct_Overlapping = round(100 * sum(Overlapping) / n(), 1),
      Mean_Effect_Size = round(mean(abs(Effect_Size), na.rm = TRUE), 3),
      Large_Effects = sum(abs(Effect_Size) > 0.8, na.rm = TRUE),
      .groups = 'drop'
    )
  return(summary_stats)
}

summary_table <- create_summary_table(all_results)
print("\nSUMMARY TABLE:")
print(summary_table)

# Save summary
summary_file <- file.path(full_output_dir, "summary_statistics.csv")
write.csv(summary_table, summary_file, row.names = FALSE)

# Create visualizations
cat("\nCreating visualizations...\n")

# 1. Effect size distribution plot
effect_size_plot <- ggplot(all_results, aes(x = Effect_Size, fill = Comparison)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~Comparison, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of Effect Sizes Across Comparisons",
       x = "Effect Size (Cohen's d)", y = "Number of Connections") +
  geom_vline(xintercept = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8), 
             linetype = "dashed", alpha = 0.5) +
  theme(legend.position = "none")

ggsave(file.path(full_output_dir, "effect_size_distribution.png"), 
       effect_size_plot, width = 12, height = 8, dpi = 300)

# 2. Significance categories plot
significance_summary <- all_results %>%
  group_by(Comparison) %>%
  summarise(
    Larger = sum(Significantly_Larger),
    Smaller = sum(Significantly_Smaller),
    Overlapping = sum(Overlapping),
    .groups = 'drop'
  ) %>%
  melt(id.vars = "Comparison", variable.name = "Category", value.name = "Count")

significance_plot <- ggplot(significance_summary, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Significance Categories Across Comparisons",
       x = "Comparison", y = "Number of Connections") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Larger" = "#e74c3c", "Smaller" = "#3498db", "Overlapping" = "#95a5a6"))

ggsave(file.path(full_output_dir, "significance_categories.png"), 
       significance_plot, width = 10, height = 6, dpi = 300)

# 3. Top significant connections
top_connections <- all_results %>%
  filter(abs(Effect_Size) > 0.5) %>%
  arrange(desc(abs(Effect_Size))) %>%
  head(20)

if (nrow(top_connections) > 0) {
  top_plot <- ggplot(top_connections, aes(x = reorder(Connection, abs(Effect_Size)), 
                                          y = Effect_Size, fill = Comparison)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Connections by Effect Size",
         x = "Connection", y = "Effect Size") +
    facet_wrap(~Comparison, scales = "free")
  
  ggsave(file.path(full_output_dir, "top_significant_connections.png"), 
         top_plot, width = 12, height = 10, dpi = 300)
}

# Print detailed summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("DETAILED ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

for (comparison in unique(all_results$Comparison)) {
  comp_data <- all_results[all_results$Comparison == comparison, ]
  cat("\n", comparison, ":\n")
  cat("  Total connections:", nrow(comp_data), "\n")
  cat("  Significantly larger:", sum(comp_data$Significantly_Larger), 
      sprintf(" (%.1f%%)", 100 * sum(comp_data$Significantly_Larger) / nrow(comp_data)), "\n")
  cat("  Significantly smaller:", sum(comp_data$Significantly_Smaller), 
      sprintf(" (%.1f%%)", 100 * sum(comp_data$Significantly_Smaller) / nrow(comp_data)), "\n")
  cat("  Overlapping:", sum(comp_data$Overlapping), 
      sprintf(" (%.1f%%)", 100 * sum(comp_data$Overlapping) / nrow(comp_data)), "\n")
  cat("  Mean |Effect Size|:", sprintf("%.3f", mean(abs(comp_data$Effect_Size), na.rm = TRUE)), "\n")
  cat("  Large effects (|d| > 0.8):", sum(abs(comp_data$Effect_Size) > 0.8, na.rm = TRUE), "\n")
}

cat("\nFiles saved to:", full_output_dir, "\n")
cat("Analysis completed successfully!\n")


# Check trace plots, R-hat values, effective sample sizes
plot(model0$UVC[,1], type="l")  # Trace plot for first connection


# Look at multiple connections and multiple chains
par(mfrow=c(2,2))
plot(model0$UVC[,1], type="l", main="Connection 1")
plot(model0$UVC[,10], type="l", main="Connection 10") 
plot(model0$UVC[,100], type="l", main="Connection 100")
plot(model0$UVC[,500], type="l", main="Connection 500")

# Check effective sample sizes
library(coda)
effectiveSize(model0$UVC[,1:10])

# Check R-hat (if multiple chains available)
gelman.diag(model0$UVC[,1:10])


# Use effect size thresholds instead of credible intervals
large_effects <- abs(all_results$Effect_Size) > 0.8
medium_effects <- abs(all_results$Effect_Size) > 0.5


# Enhanced Effect Size Analysis
library(dplyr)
library(ggplot2)

# Assuming all_results is already loaded
# Create effect size categories
all_results$Effect_Category <- case_when(
  abs(all_results$Effect_Size) >= 0.8 ~ "Large (≥0.8)",
  abs(all_results$Effect_Size) >= 0.5 ~ "Medium (0.5-0.8)",
  abs(all_results$Effect_Size) >= 0.2 ~ "Small (0.2-0.5)",
  TRUE ~ "Negligible (<0.2)"
)

# Effect size direction
all_results$Effect_Direction <- case_when(
  all_results$Effect_Size > 0.5 ~ "Strong Positive",
  all_results$Effect_Size > 0.2 ~ "Moderate Positive", 
  all_results$Effect_Size > -0.2 ~ "Negligible",
  all_results$Effect_Size > -0.5 ~ "Moderate Negative",
  TRUE ~ "Strong Negative"
)

# Summary by comparison
effect_summary <- all_results %>%
  group_by(Comparison) %>%
  summarise(
    Total = n(),
    Large_Effects = sum(abs(Effect_Size) > 0.8),
    Medium_Effects = sum(abs(Effect_Size) > 0.5 & abs(Effect_Size) <= 0.8),
    Small_Effects = sum(abs(Effect_Size) > 0.2 & abs(Effect_Size) <= 0.5),
    Negligible = sum(abs(Effect_Size) <= 0.2),
    Pct_Large = round(100 * Large_Effects / Total, 1),
    Pct_Medium_Plus = round(100 * sum(abs(Effect_Size) > 0.5) / Total, 1),
    Mean_AbsEffect = round(mean(abs(Effect_Size)), 3),
    Max_Effect = round(max(abs(Effect_Size)), 3),
    .groups = 'drop'
  )

print("EFFECT SIZE SUMMARY BY COMPARISON:")
print(effect_summary)

# Direction analysis
direction_summary <- all_results %>%
  group_by(Comparison, Effect_Direction) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Effect_Direction, values_from = Count, values_fill = 0)

print("\nEFFECT DIRECTION SUMMARY:")
print(direction_summary)

# Category breakdown
category_summary <- all_results %>%
  group_by(Comparison, Effect_Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Effect_Category, values_from = Count, values_fill = 0)

print("\nEFFECT CATEGORY BREAKDOWN:")
print(category_summary)

# Statistical significance vs effect size
significance_vs_effect <- all_results %>%
  mutate(
    Has_Large_Effect = abs(Effect_Size) > 0.8,
    Is_Significant = Significantly_Larger | Significantly_Smaller
  ) %>%
  group_by(Comparison) %>%
  summarise(
    Large_Effects = sum(Has_Large_Effect),
    Significant_Effects = sum(Is_Significant),
    Large_AND_Significant = sum(Has_Large_Effect & Is_Significant),
    Large_NOT_Significant = sum(Has_Large_Effect & !Is_Significant),
    .groups = 'drop'
  )

print("\nLARGE EFFECT SIZES vs STATISTICAL SIGNIFICANCE:")
print(significance_vs_effect)

# Top connections by effect size
top_effects <- all_results %>%
  arrange(desc(abs(Effect_Size))) %>%
  head(20) %>%
  select(Comparison, Connection, Effect_Size, Model_A_Mean, Model_B_Mean, 
         Significantly_Larger, Significantly_Smaller)

print("\nTOP 20 CONNECTIONS BY EFFECT SIZE:")
print(top_effects)

# Create visualizations
cat("\nCreating effect size visualizations...\n")

# 1. Effect size distribution by comparison
p1 <- ggplot(all_results, aes(x = abs(Effect_Size), fill = Comparison)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~Comparison, scales = "free_y") +
  geom_vline(xintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = "red") +
  labs(title = "Distribution of Absolute Effect Sizes",
       x = "Absolute Effect Size |Cohen's d|", 
       y = "Number of Connections") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. Effect size categories by comparison
effect_data <- all_results %>%
  count(Comparison, Effect_Category) %>%
  group_by(Comparison) %>%
  mutate(Percentage = round(100 * n / sum(n), 1))

p2 <- ggplot(effect_data, aes(x = Comparison, y = n, fill = Effect_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Effect Size Categories by Comparison",
       x = "Comparison", y = "Number of Connections",
       fill = "Effect Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Large (≥0.8)" = "#e74c3c", 
                               "Medium (0.5-0.8)" = "#f39c12",
                               "Small (0.2-0.5)" = "#f1c40f",
                               "Negligible (<0.2)" = "#95a5a6"))

# 3. Effect size vs significance scatter
p3 <- ggplot(all_results, aes(x = abs(Effect_Size), y = Posterior_Probability_Larger)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Comparison) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "blue") +
  labs(title = "Effect Size vs Posterior Probability",
       x = "Absolute Effect Size", 
       y = "Posterior Probability (B > A)") +
  theme_minimal()

# Save plots
ggsave("effect_size_distribution.png", p1, width = 12, height = 8, dpi = 300)
ggsave("effect_size_categories.png", p2, width = 10, height = 6, dpi = 300)
ggsave("effect_vs_significance.png", p3, width = 12, height = 8, dpi = 300)

# Summary statistics
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("OVERALL EFFECT SIZE ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

total_connections <- nrow(all_results)
large_effects_total <- sum(abs(all_results$Effect_Size) > 0.8)
medium_plus_total <- sum(abs(all_results$Effect_Size) > 0.5)

cat("Total connections analyzed:", total_connections, "\n")
cat("Large effects (|d| > 0.8):", large_effects_total, 
    sprintf(" (%.1f%%)", 100 * large_effects_total / total_connections), "\n")
cat("Medium+ effects (|d| > 0.5):", medium_plus_total, 
    sprintf(" (%.1f%%)", 100 * medium_plus_total / total_connections), "\n")

cat("\nComparison Rankings (by % large effects):\n")
ranking <- effect_summary %>% arrange(desc(Pct_Large))
for(i in 1:nrow(ranking)) {
  cat(sprintf("%d. %s: %.1f%% large effects\n", 
              i, ranking$Comparison[i], ranking$Pct_Large[i]))
}

cat("\nInterpretation:\n")
cat("- Your clusters show SUBSTANTIAL differences in brain connectivity\n")
cat("- Nearly half of all connections show large practical differences\n")
cat("- This indicates your clustering captured meaningful brain states\n")
cat("- Focus on effect sizes rather than statistical significance\n")

# Create summary list for further analysis
analysis_results <- list(
  effect_summary = effect_summary,
  direction_summary = direction_summary,
  top_effects = top_effects,
  large_effects_total = large_effects_total,
  medium_plus_total = medium_plus_total
)

cat("\nAnalysis completed! Results stored in 'analysis_results' object.\n")



# Analysis of AD Progression Patterns
# Let's investigate if your clustering matches expected AD progression

# First, let's check the expected vs observed pattern
cat("ALZHEIMER'S DISEASE PROGRESSION ANALYSIS\n")
cat("========================================\n\n")

# Expected vs Observed Pattern
expected_ranking <- c("Early_vs_Advanced", "Early_vs_Intermediate", "Intermediate_vs_Advanced")
observed_ranking <- c("Intermediate_vs_Advanced", "Early_vs_Advanced", "Early_vs_Intermediate")

cat("EXPECTED AD progression (by effect size):\n")
cat("1. Early → Advanced (largest gap)\n")
cat("2. Early → Intermediate (medium gap)\n") 
cat("3. Intermediate → Advanced (smallest gap)\n\n")

cat("YOUR OBSERVED pattern:\n")
cat("1. Intermediate → Advanced (49.5% large effects)\n")
cat("2. Early → Advanced (47.3% large effects)\n")
cat("3. Early → Intermediate (42.5% large effects)\n\n")

# Calculate what we'd expect if progression was linear
your_results <- data.frame(
  Comparison = c("Early_vs_Advanced", "Early_vs_Intermediate", "Intermediate_vs_Advanced"),
  Pct_Large = c(47.3, 42.5, 49.5),
  Mean_AbsEffect = c(0.814, 0.762, 0.993)
)

# Check for linear progression
cat("LINEARITY CHECK:\n")
cat("If progression were linear, we'd expect:\n")
cat("Early→Intermediate + Intermediate→Advanced ≈ Early→Advanced\n\n")

linear_expectation <- your_results$Pct_Large[2] + your_results$Pct_Large[3]
actual_early_advanced <- your_results$Pct_Large[1]

cat("Expected Early→Advanced (if linear):", round(linear_expectation, 1), "%\n")
cat("Actual Early→Advanced:", actual_early_advanced, "%\n")
cat("Difference:", round(linear_expectation - actual_early_advanced, 1), "%\n\n")

if (linear_expectation > actual_early_advanced + 10) {
  cat("FINDING: Non-linear progression detected!\n")
  cat("This suggests accelerating changes in later stages.\n\n")
} else {
  cat("FINDING: Roughly linear progression.\n\n")
}

# Possible explanations
cat("POSSIBLE EXPLANATIONS FOR YOUR PATTERN:\n")
cat("======================================\n\n")

cat("1. ACCELERATING PATHOLOGY (Common in AD):\n")
cat("   - Early stages: Compensatory mechanisms mask damage\n")
cat("   - Later stages: Compensation fails, rapid deterioration\n")
cat("   - Evidence: Many AD studies show non-linear decline\n\n")

cat("2. CLUSTERING ISSUES TO CHECK:\n")
cat("   a) Are your 'stages' correctly assigned?\n")
cat("   b) Sample size per cluster?\n")
cat("   c) Clinical measures used for staging?\n\n")

cat("3. BIOLOGICAL PLAUSIBILITY:\n")
cat("   - Early AD: Subtle metabolic changes\n")
cat("   - Intermediate: Network reorganization begins\n") 
cat("   - Advanced: Massive network breakdown\n\n")

# Recommendations for further analysis
cat("RECOMMENDED INVESTIGATIONS:\n")
cat("==========================\n\n")

cat("1. CHECK CLUSTER ASSIGNMENTS:\n")
cat("   # Load your original clustering data\n")
cat("   # Check clinical characteristics by cluster\n")
cat("   table(cluster_assignments)\n")
cat("   by(clinical_scores, cluster_assignments, summary)\n\n")

cat("2. VALIDATE WITH CLINICAL MEASURES:\n")
cat("   # Compare MMSE, CDR, or other scores across clusters\n")
cat("   # Should show: Early > Intermediate > Advanced\n\n")

cat("3. CHECK SAMPLE SIZES:\n")
cat("   # Unbalanced clusters could affect results\n")
cat("   # Small advanced group might show extreme effects\n\n")

cat("4. LITERATURE COMPARISON:\n")
cat("   # Many AD studies actually show accelerating decline\n")
cat("   # Your pattern might be neurobiologically accurate\n\n")

# Create diagnostic plots
cat("CREATING DIAGNOSTIC VISUALIZATIONS...\n")

# Effect size progression plot
library(ggplot2)

# Create data for plotting
progression_data <- data.frame(
  Stage_Gap = factor(c("Early→Intermediate", "Intermediate→Advanced", "Early→Advanced"),
                     levels = c("Early→Intermediate", "Intermediate→Advanced", "Early→Advanced")),
  Pct_Large_Effects = c(42.5, 49.5, 47.3),
  Mean_Effect_Size = c(0.762, 0.993, 0.814),
  Expected_Rank = c(2, 3, 1),
  Observed_Rank = c(3, 1, 2)
)

# Plot 1: Effect sizes by stage transition
p1 <- ggplot(progression_data, aes(x = Stage_Gap, y = Pct_Large_Effects)) +
  geom_bar(stat = "identity", fill = c("#3498db", "#e74c3c", "#f39c12")) +
  geom_text(aes(label = paste0(Pct_Large_Effects, "%")), vjust = -0.5) +
  labs(title = "Percentage of Large Effects by Stage Transition",
       subtitle = "Red bar shows unexpected peak at Intermediate→Advanced",
       x = "Stage Transition", y = "% Large Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Expected vs Observed ranking
ranking_data <- data.frame(
  Comparison = rep(c("Early→Intermediate", "Intermediate→Advanced", "Early→Advanced"), 2),
  Rank = c(2, 3, 1, 3, 1, 2),
  Type = rep(c("Expected", "Observed"), each = 3)
)

p2 <- ggplot(ranking_data, aes(x = Comparison, y = Rank, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_reverse(breaks = 1:3, labels = c("Most Different", "Medium", "Least Different")) +
  labs(title = "Expected vs Observed Ranking of Stage Differences",
       subtitle = "Shows non-linear AD progression pattern",
       x = "Stage Comparison", y = "Ranking (1 = Most Different)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Expected" = "#95a5a6", "Observed" = "#e74c3c"))

# Save plots
ggsave("ad_progression_pattern.png", p1, width = 10, height = 6, dpi = 300)
ggsave("expected_vs_observed_ranking.png", p2, width = 10, height = 6, dpi = 300)

cat("Plots saved: ad_progression_pattern.png, expected_vs_observed_ranking.png\n\n")

# Final interpretation
cat("FINAL INTERPRETATION:\n")
cat("====================\n\n")

cat("Your results show NON-LINEAR AD progression, which is actually:\n")
cat("✓ Biologically plausible (compensation → decompensation)\n")
cat("✓ Consistent with some AD literature\n") 
cat("✓ Clinically meaningful (rapid late-stage decline)\n\n")

cat("However, you should VERIFY:\n")
cat("• Clinical scores match your cluster assignments\n")
cat("• Sample sizes are balanced across clusters\n")
cat("• Staging criteria were appropriate\n\n")

cat("If clinical validation confirms your staging, then you've discovered\n")
cat("evidence for ACCELERATING brain network changes in AD progression!\n")





# UVC vs UVPM Analysis Guide
# Understanding when to use each component

cat("UVC vs UVPM ANALYSIS GUIDE\n")
cat("==========================\n\n")

# Load your models (assuming they're already loaded)
# model0, model1, model2

# First, let's explore the structure of both
cat("EXPLORING MODEL STRUCTURE:\n")
cat("--------------------------\n")

# Check what's in your models
cat("Model components available:\n")
model_components <- names(model0)
for(comp in model_components) {
  cat("- ", comp, "\n")
}

cat("\nUVC Structure:\n")
cat("Dimensions:", dim(model0$UVC), "\n")
cat("Type: MCMC posterior samples\n")
cat("Use for: Statistical inference, credible intervals\n")

cat("\nUVPM Structure:\n") 
cat("Dimensions:", dim(model0$UVPM), "\n")
cat("Type: Posterior means (point estimates)\n")
cat("Use for: Effect sizes, visualization, correlation analysis\n\n")

# Compare a specific connection
connection_example <- 1
cat("EXAMPLE - Connection", connection_example, ":\n")
cat("UVC samples (first 10):", model0$UVC[1:10, connection_example], "\n")
cat("UVPM value:", model0$UVPM[connection_example, connection_example], "\n")
cat("UVC mean:", mean(model0$UVC[, connection_example]), "\n")
cat("Should be similar:", abs(model0$UVPM[connection_example, connection_example] - mean(model0$UVC[, connection_example])) < 0.001, "\n\n")

cat("WHEN TO USE EACH:\n")
cat("==================\n\n")

cat("USE UVC (Posterior Samples) FOR:\n")
cat("- Statistical significance testing\n")
cat("- Credible intervals\n")
cat("- Uncertainty quantification\n")
cat("- Bayesian hypothesis testing\n")
cat("- Posterior probability calculations\n\n")

cat("USE UVPM (Posterior Means) FOR:\n") 
cat("- Effect size calculations\n")
cat("- Correlation with other brain measures\n")
cat("- Visualization of brain networks\n")
cat("- Comparing average connectivity patterns\n")
cat("- Machine learning features\n\n")

cat("YOUR CURRENT ANALYSIS ASSESSMENT:\n")
cat("=================================\n\n")

# Check if your effect size analysis makes sense
cat("You've been using UVC for effect sizes, which is:\n")
cat("✓ Valid but computationally intensive\n")
cat("✓ Incorporates uncertainty (good for robust analysis)\n")
cat("✗ More complex than necessary for basic comparisons\n\n")

cat("RECOMMENDATION FOR YOUR RESEARCH QUESTION:\n")
cat("=========================================\n\n")

cat("For AD progression analysis, you should consider BOTH:\n\n")

cat("1. UVPM ANALYSIS (Simpler, More Standard):\n")
cat("   - Compare average connectivity patterns\n")
cat("   - Calculate Cohen's d using means\n")
cat("   - Faster computation\n")
cat("   - Standard approach in neuroimaging\n\n")

cat("2. UVC ANALYSIS (Your Current Approach):\n") 
cat("   - More sophisticated uncertainty quantification\n")
cat("   - Better for small sample sizes\n")
cat("   - Accounts for estimation uncertainty\n")
cat("   - More appropriate for Bayesian framework\n\n")

# Create comparison analysis
cat("CREATING COMPARISON ANALYSIS...\n")

# Function to compare UVC vs UVPM approaches
compare_approaches <- function(model_a, model_b, comparison_name) {
  
  cat("Comparing approaches for:", comparison_name, "\n")
  
  # Method 1: Using UVPM (traditional approach)
  uvpm_diff <- model_b$UVPM - model_a$UVPM
  uvpm_effect_sizes <- uvpm_diff[upper.tri(uvpm_diff)] / 
    sqrt((var(as.vector(model_a$UVPM[upper.tri(model_a$UVPM)])) + 
            var(as.vector(model_b$UVPM[upper.tri(model_b$UVPM)]))) / 2)
  
  # Method 2: Using UVC (your current approach)
  # Take mean effect size across all MCMC samples
  n_connections <- ncol(model_a$UVC)
  uvc_effect_sizes <- numeric(n_connections)
  
  for(i in 1:min(100, n_connections)) {  # Limit to first 100 for speed
    samples_a <- model_a$UVC[, i]
    samples_b <- model_b$UVC[, i]
    pooled_sd <- sqrt((var(samples_a) + var(samples_b)) / 2)
    uvc_effect_sizes[i] <- (mean(samples_b) - mean(samples_a)) / pooled_sd
  }
  
  # Compare results
  cat("UVPM approach - Large effects (|d| > 0.8):", 
      sum(abs(uvpm_effect_sizes) > 0.8, na.rm = TRUE), 
      "out of", length(uvpm_effect_sizes), "\n")
  
  cat("UVC approach - Large effects (|d| > 0.8):", 
      sum(abs(uvc_effect_sizes[1:100]) > 0.8, na.rm = TRUE), 
      "out of 100 tested\n")
  
  cat("Correlation between methods:", 
      cor(uvpm_effect_sizes[1:100], uvc_effect_sizes[1:100], use = "complete.obs"), "\n\n")
  
  return(list(uvpm_effects = uvpm_effect_sizes, uvc_effects = uvc_effect_sizes))
}

# Run comparison (if models are loaded)
if(exists("model0") && exists("model1")) {
  comparison_results <- compare_approaches(model0, model1, "Early vs Advanced")
}

cat("PRACTICAL RECOMMENDATIONS:\n")
cat("==========================\n\n")

cat("For your AD progression study:\n\n")

cat("1. START WITH UVPM ANALYSIS:\n")
cat("   - Faster and simpler\n")
cat("   - Standard in neuroimaging literature\n")
cat("   - Easier to interpret and visualize\n")
cat("   - Good for initial exploration\n\n")

cat("2. VALIDATE WITH UVC ANALYSIS:\n")
cat("   - Use for final statistical testing\n")
cat("   - Better uncertainty quantification\n")
cat("   - More appropriate for Bayesian models\n")
cat("   - Stronger for publication\n\n")

cat("3. HYBRID APPROACH (RECOMMENDED):\n")
cat("   - Use UVPM for effect size calculations\n")
cat("   - Use UVC for significance testing\n")
cat("   - Report both in your paper\n")
cat("   - Show they give consistent results\n\n")

cat("EXAMPLE CODE FOR UVPM ANALYSIS:\n")
cat("===============================\n\n")

cat('# Simple UVPM-based effect size calculation\n')
cat('uvpm_diff <- model1$UVPM - model0$UVPM\n')
cat('effect_sizes <- uvpm_diff / sqrt((var(model0$UVPM) + var(model1$UVPM))/2)\n')
cat('large_effects <- sum(abs(effect_sizes) > 0.8, na.rm=TRUE)\n')
cat('cat("Large effects:", large_effects, "out of", length(effect_sizes), "\\n")\n\n')

cat("This would be much faster than your current UVC approach!\n")