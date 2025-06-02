# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)
library(RColorBrewer)
theme_minimal()

# Read your data
pet_data <- read.csv("C:/Users/theya/Downloads/ADNIFinal/Brainplotk=3.csv")

# 1. OVERLAID GROUP HISTOGRAMS/DENSITY CURVES
# Transform data to long format for plotting
pet_long <- pet_data %>%
  select(PredictedLabel3, starts_with("V")) %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "region", 
               values_to = "tau_burden") %>%
  mutate(cluster = factor(paste("Cluster", PredictedLabel3)))

# Plot 1: Overlaid histograms with density curves for all regions combined
p1 <- ggplot(pet_long, aes(x = tau_burden, fill = cluster, color = cluster)) +
  geom_histogram(aes(y = after_stat(density)), 
                 alpha = 0.6, bins = 30, position = "identity") +
  geom_density(alpha = 0.8, size = 1.2) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  labs(title = "Distribution of Regional Tau Burden Across Clusters",
       subtitle = "All 68 brain regions combined",
       x = "Tau Burden (SUVR)",
       y = "Density",
       fill = "Cluster",
       color = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "top")

# Plot 2: Separate density plots for key regions (you can modify regions)
key_regions <- c("V26", "V27", "V60", "V61")  # parahippocampal and entorhinal (both hemispheres)
region_names <- c("Parahippocampal R", "Entorhinal R", "Parahippocampal L", "Entorhinal L")

pet_key_regions <- pet_data %>%
  select(PredictedLabel3, all_of(key_regions)) %>%
  pivot_longer(cols = all_of(key_regions), 
               names_to = "region", 
               values_to = "tau_burden") %>%
  mutate(cluster = factor(paste("Cluster", PredictedLabel3)),
         region_label = factor(region, 
                               levels = key_regions, 
                               labels = region_names))

p2 <- ggplot(pet_key_regions, aes(x = tau_burden, fill = cluster)) +
  geom_histogram(aes(y = after_stat(density)), 
                 alpha = 0.7, bins = 20, position = "identity") +
  geom_density(aes(color = cluster), alpha = 0, size = 1) +
  facet_wrap(~region_label, scales = "free", nrow = 2) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  labs(title = "Regional Tau Distribution in Key AD-Related Regions",
       x = "Tau Burden (SUVR)",
       y = "Density") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

# 2. INDIVIDUAL SUBJECT DENSITY CURVES
# Create individual subject curves for a subset of subjects (to avoid overcrowding)
set.seed(123)
sample_subjects <- pet_data %>%
  group_by(PredictedLabel3) %>%
  slice_sample(n = 10) %>%  # Sample 10 subjects per cluster
  ungroup()

# Transform sample data
sample_long <- sample_subjects %>%
  select(PredictedLabel3, starts_with("V")) %>%
  mutate(subject_id = row_number()) %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "region", 
               values_to = "tau_burden") %>%
  mutate(cluster = factor(paste("Cluster", PredictedLabel3)))

# Plot 3: Individual subject density curves
p3 <- ggplot() +
  # Individual subject curves (lighter)
  stat_density(data = sample_long, 
               aes(x = tau_burden, group = subject_id, color = cluster),
               alpha = 0.3, size = 0.5, geom = "line") +
  # Group mean curves (bolder)
  stat_density(data = pet_long, 
               aes(x = tau_burden, color = cluster),
               alpha = 1, size = 2, geom = "line") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  labs(title = "Individual Subject Density Curves by Cluster",
       subtitle = "Light lines: individual subjects, Bold lines: group means",
       x = "Tau Burden (SUVR)",
       y = "Density",
       color = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"))

# 3. COMPREHENSIVE SUMMARY STATISTICS
summary_stats <- pet_long %>%
  group_by(cluster) %>%
  summarise(
    n_observations = n(),
    mean_tau = mean(tau_burden, na.rm = TRUE),
    median_tau = median(tau_burden, na.rm = TRUE),
    sd_tau = sd(tau_burden, na.rm = TRUE),
    q25 = quantile(tau_burden, 0.25, na.rm = TRUE),
    q75 = quantile(tau_burden, 0.75, na.rm = TRUE),
    min_tau = min(tau_burden, na.rm = TRUE),
    max_tau = max(tau_burden, na.rm = TRUE),
    .groups = "drop"
  )

print("Summary Statistics by Cluster:")
print(summary_stats)

# 4. VIOLIN PLOT FOR COMPARISON
p4 <- ggplot(pet_long, aes(x = cluster, y = tau_burden, fill = cluster)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "Tau Burden Distribution by Cluster",
       subtitle = "Violin plots showing full distribution shape",
       x = "Cluster",
       y = "Tau Burden (SUVR)") +
  theme_minimal() +
  theme(legend.position = "none")

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)

# Save plots
ggsave(file.path(output_dir, "tau_distributions_overlay.png"), p1, 
       width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "tau_distributions_key_regions.png"), p2, 
       width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "individual_density_curves.png"), p3, 
       width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "tau_violin_plots.png"), p4, 
       width = 10, height = 8, dpi = 300)

# 5. ADVANCED: Regional comparison across clusters
# Calculate mean tau burden per region per cluster
regional_summary <- pet_data %>%
  select(PredictedLabel3, starts_with("V")) %>%
  group_by(PredictedLabel3) %>%
  summarise(across(starts_with("V"), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "region", 
               values_to = "mean_tau") %>%
  mutate(cluster = factor(paste("Cluster", PredictedLabel3)))

# Heatmap of regional differences
p5 <- ggplot(regional_summary, aes(x = region, y = cluster, fill = mean_tau)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Tau\nBurden") +
  labs(title = "Regional Tau Burden Heatmap Across Clusters",
       x = "Brain Region",
       y = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

print(p5)
ggsave(file.path(output_dir, "regional_tau_heatmap.png"), p5, 
       width = 16, height = 6, dpi = 300)

# Statistical tests between clusters
# Perform ANOVA for overall differences
aov_result <- aov(tau_burden ~ cluster, data = pet_long)
print("ANOVA Results:")
print(summary(aov_result))

# Pairwise t-tests
pairwise_results <- pairwise.t.test(pet_long$tau_burden, pet_long$cluster, 
                                    p.adjust.method = "bonferroni")
print("Pairwise t-test results:")
print(pairwise_results)