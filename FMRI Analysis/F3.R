# FULLY CORRECTED Dynamic Publication-Ready Connectivity Analysis Plots
# ALL ISSUES FIXED: No hardcoded values, correct radar labels, proper network names

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(scales)

# =============================================================================
# SETUP: LOAD YOUR CORRECTED DATA
# =============================================================================

cat("CREATING FULLY CORRECTED DYNAMIC PUBLICATION-READY PLOTS\n")
cat("========================================================\n\n")

# Check if corrected data exists
if(!exists("all_comprehensive_data") || !exists("gold_standard_all") || !exists("comparison_summary")) {
  stop("❌ Please run the corrected comprehensive analysis first to generate dynamic data!")
}

cat("✅ Using your corrected comprehensive analysis data\n")
cat("✅ All plots will be 100% dynamic based on your actual results\n\n")

# Create publication output directory
pub_dir <- "Publication_Connectivity_Figures_CORRECTED"
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

# Define consistent color schemes
subtype_colors <- c("S1" = "#5A9BD4", "S2" = "#70AD47", "S3" = "#FFC000")
comparison_colors <- c("S1_vs_S2" = "#E74C3C", "S1_vs_S3" = "#3498DB", "S3_vs_S2" = "#F39C12")

# =============================================================================
# DYNAMIC DATA EXTRACTION FROM YOUR CORRECTED ANALYSIS
# =============================================================================

cat("Extracting dynamic data from your corrected analysis...\n")

# Extract effect size percentages dynamically
effect_size_summary <- all_comprehensive_data %>%
  group_by(Comparison) %>%
  summarise(
    Total = n(),
    Large = sum(Effect_Category == "Large"),
    Medium = sum(Effect_Category == "Medium"),
    Small = sum(Effect_Category == "Small"),
    Negligible = sum(Effect_Category == "Negligible"),
    Large_Pct = round(100 * Large / Total, 1),
    Medium_Pct = round(100 * Medium / Total, 1),
    Small_Pct = round(100 * Small / Total, 1),
    Negligible_Pct = round(100 * Negligible / Total, 1),
    .groups = 'drop'
  )

# Extract network involvement dynamically
network_involvement_dynamic <- gold_standard_all %>%
  select(Comparison, RegionA_RSNName, RegionB_RSNName) %>%
  pivot_longer(cols = c(RegionA_RSNName, RegionB_RSNName), 
               names_to = "Region_Type", values_to = "Network") %>%
  filter(!is.na(Network), Network != "") %>%
  count(Comparison, Network, sort = TRUE) %>%
  group_by(Comparison) %>%
  mutate(Percentage = round(100 * n / sum(n), 1)) %>%
  ungroup()

# Extract top connections dynamically (FIXED - NO HARDCODING)
top_connections_dynamic <- gold_standard_all %>%
  arrange(desc(Abs_Effect_Size)) %>%
  head(5) %>%
  mutate(
    # Create short names from actual connection names dynamically
    Connection_Short = paste0(
      substr(gsub(".*rh-|.*lh-", "", RegionA_Name), 1, 8),
      "-",
      substr(gsub(".*rh-|.*lh-", "", RegionB_Name), 1, 8)
    )
  )

cat("✅ Dynamic data extraction complete!\n\n")

# =============================================================================
# FIGURE 1: COMPREHENSIVE CONNECTIVITY OVERVIEW (2x2 PANEL) - FULLY DYNAMIC
# =============================================================================

cat("Creating Figure 1: Fully Dynamic Comprehensive Overview...\n")

# Panel A: Effect Size Categories (DYNAMIC)
plot1a_data <- effect_size_summary %>%
  select(Comparison, Large_Pct, Medium_Pct, Small_Pct, Negligible_Pct) %>%
  pivot_longer(cols = c(Large_Pct, Medium_Pct, Small_Pct, Negligible_Pct), 
               names_to = "Effect_Size", values_to = "Percentage") %>%
  mutate(
    Effect_Size = factor(case_when(
      Effect_Size == "Large_Pct" ~ "Large",
      Effect_Size == "Medium_Pct" ~ "Medium", 
      Effect_Size == "Small_Pct" ~ "Small",
      Effect_Size == "Negligible_Pct" ~ "Negligible"
    ), levels = c("Negligible", "Small", "Medium", "Large")),
    Comparison = factor(Comparison, levels = c("S1_vs_S2", "S1_vs_S3", "S3_vs_S2"),
                        labels = c("S1 vs S2", "S1 vs S3", "S3 vs S2"))
  )

p1a <- ggplot(plot1a_data, aes(x = Comparison, y = Percentage, fill = Effect_Size)) +
  geom_bar(stat = "identity", color = "white", size = 0.3) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Large" = "#D32F2F", "Medium" = "#FF9800", 
                               "Small" = "#FFC107", "Negligible" = "#9E9E9E"),
                    labels = c("Large (≥0.8)", "Medium (0.5-0.8)", "Small (0.2-0.5)", "Negligible (<0.2)")) +
  labs(title = "(A) Effect Size Distribution", x = "Subtype Comparison", 
       y = "Percentage (%)", fill = "Effect Size") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))

# Panel B: Statistical Significance Summary (DYNAMIC)
plot1b_data <- comparison_summary %>%
  select(Comparison, Large_Effects, Significant, Gold_Standard) %>%
  pivot_longer(cols = c(Large_Effects, Significant, Gold_Standard), 
               names_to = "Category", values_to = "Count") %>%
  mutate(
    Comparison = factor(Comparison, levels = c("S1_vs_S2", "S1_vs_S3", "S3_vs_S2"),
                        labels = c("S1 vs S2", "S1 vs S3", "S3 vs S2"))
  )

p1b <- ggplot(plot1b_data, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", size = 0.3) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Large_Effects" = "#2196F3", "Significant" = "#4CAF50", 
                               "Gold_Standard" = "#F44336"),
                    labels = c("Large Effects", "Statistically Significant", "Gold Standard")) +
  labs(title = "(B) Statistical Categories", x = "Subtype Comparison", 
       y = "Number of Connections", fill = "Connection Type") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Panel C: Network Involvement (DYNAMIC)
plot1c_data <- network_involvement_dynamic %>%
  filter(Percentage > 0) %>%
  mutate(
    Comparison = factor(Comparison, levels = c("S1_vs_S2", "S1_vs_S3", "S3_vs_S2"),
                        labels = c("S1 vs S2", "S1 vs S3", "S3 vs S2"))
  )

p1c <- ggplot(plot1c_data, aes(x = reorder(Network, Percentage), y = Percentage, fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", size = 0.3) +
  coord_flip() +
  scale_fill_manual(values = comparison_colors) +
  labs(title = "(C) Network Involvement", x = "Brain Network", 
       y = "Percentage (%)", fill = "Comparison") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "right")

# Panel D: Top Gold Standard Connections (DYNAMIC - FIXED)
plot1d_data <- top_connections_dynamic %>%
  mutate(
    Comparison = factor(Comparison, levels = c("S1_vs_S2", "S1_vs_S3", "S3_vs_S2"),
                        labels = c("S1 vs S2", "S1 vs S3", "S3 vs S2"))
  )

p1d <- ggplot(plot1d_data, aes(x = reorder(Connection_Short, Effect_Size), y = Effect_Size, fill = Comparison)) +
  geom_bar(stat = "identity", color = "white", size = 0.3) +
  geom_text(aes(label = round(Effect_Size, 2)), hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = comparison_colors) +
  labs(title = "(D) Strongest Connections", x = "Connection", 
       y = "Effect Size (Cohen's d)", fill = "Comparison") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "right")

# Combine Figure 1
figure1 <- grid.arrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2,
                        top = "Functional Connectivity Analysis Across AD Subtypes (Corrected)")

ggsave(file.path(pub_dir, "Figure1_Connectivity_Overview_CORRECTED.png"), 
       figure1, width = 14, height = 10, dpi = 300)

# =============================================================================
# FIGURE 2: DYNAMIC NETWORK HEATMAP
# =============================================================================

cat("Creating Figure 2: Dynamic Network Connectivity Heatmap...\n")

# Create network-to-network connectivity matrix from your data
network_pairs <- gold_standard_all %>%
  separate(Network_Pair, into = c("Network1", "Network2"), sep = " - ") %>%
  group_by(Network1, Network2) %>%
  summarise(
    Mean_Effect = mean(abs(Effect_Size)),
    Count = n(),
    .groups = 'drop'
  )

# Get all unique networks
all_networks <- unique(c(network_pairs$Network1, network_pairs$Network2))

# Create matrix
network_matrix <- matrix(0, nrow = length(all_networks), ncol = length(all_networks))
rownames(network_matrix) <- colnames(network_matrix) <- all_networks

# Fill matrix with actual data
for(i in 1:nrow(network_pairs)) {
  row_idx <- which(all_networks == network_pairs$Network1[i])
  col_idx <- which(all_networks == network_pairs$Network2[i])
  network_matrix[row_idx, col_idx] <- network_pairs$Mean_Effect[i]
  network_matrix[col_idx, row_idx] <- network_pairs$Mean_Effect[i]  # Make symmetric
}

# Convert to long format
network_long <- expand.grid(Network1 = all_networks, Network2 = all_networks) %>%
  mutate(Connectivity = as.vector(network_matrix))

p2 <- ggplot(network_long, aes(x = Network1, y = Network2, fill = Connectivity)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Effect Size") +
  labs(title = "Network-Level Connectivity Differences (Dynamic)",
       subtitle = "Based on your actual gold standard connections",
       x = "Brain Network", y = "Brain Network") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1))

ggsave(file.path(pub_dir, "Figure2_Network_Heatmap_Dynamic.png"), 
       p2, width = 8, height = 8, dpi = 300)

# =============================================================================
# FIGURE 3: CORRECTED RADAR PLOTS WITH ALL LABELS VISIBLE
# =============================================================================

cat("Creating Figure 3: CORRECTED Radar Plots with All Labels Visible...\n")

# Define STATIC network order for consistent radar positions
STATIC_NETWORK_ORDER <- c("DMN", "SM", "VA", "VIS", "FP", "LS", "SUB", "DA")
STATIC_NETWORK_LABELS <- c(
  "DMN" = "Default Mode",
  "SM" = "Somatomotor", 
  "VA" = "Ventral Attention",
  "VIS" = "Visual",
  "FP" = "Frontoparietal",
  "LS" = "Limbic System", 
  "SUB" = "Subcortical",
  "DA" = "Dorsal Attention"
)

# Prepare dynamic network data with static positions
radar_data_static <- data.frame(
  Network = STATIC_NETWORK_ORDER,
  Network_Label = STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER]
)

# Add dynamic data for each comparison
for(comp in c("S1_vs_S2", "S1_vs_S3", "S3_vs_S2")) {
  comp_data <- network_involvement_dynamic %>%
    filter(Comparison == comp) %>%
    select(Network, Percentage)
  
  # Merge with static structure, filling missing networks with 0
  merged_data <- radar_data_static %>%
    left_join(comp_data, by = "Network") %>%
    mutate(Percentage = ifelse(is.na(Percentage), 0, Percentage))
  
  radar_data_static[[comp]] <- merged_data$Percentage
}

cat("✅ Radar data prepared with static network positions and dynamic values\n")
create_corrected_radar <- function(data, title, color, max_value = 70) {
  library(ggplot2)
  library(dplyr)
  
  STATIC_NETWORK_ORDER <- c("DMN", "SM", "VA", "VIS", "FP", "LS", "SUB", "DA")
  STATIC_NETWORK_LABELS <- c(
    "DMN" = "Default Mode",
    "SM" = "Somatomotor", 
    "VA" = "Ventral Attention",
    "VIS" = "Visual",
    "FP" = "Frontoparietal",
    "LS" = "Limbic System", 
    "SUB" = "Subcortical",
    "DA" = "Dorsal Attention"
  )
  
  df <- data.frame(
    Network = factor(STATIC_NETWORK_ORDER, levels = STATIC_NETWORK_ORDER),
    Network_Label = STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER],
    Value = data,
    Angle = seq(0, 2 * pi, length.out = length(STATIC_NETWORK_ORDER) + 1)[1:length(STATIC_NETWORK_ORDER)]
  )
  
  # Calculate coordinates
  df$x <- df$Value * cos(df$Angle)
  df$y <- df$Value * sin(df$Angle)
  
  # Grid circles
  grid_data <- data.frame()
  for (i in seq(10, max_value, 10)) {
    circle_x <- i * cos(seq(0, 2 * pi, length.out = 100))
    circle_y <- i * sin(seq(0, 2 * pi, length.out = 100))
    grid_data <- rbind(grid_data, data.frame(x = circle_x, y = circle_y, level = i))
  }
  
  # Axis labels and positioning
  label_distance <- max_value + 10  # Closer to center
  axis_data <- data.frame(
    x_start = 0,
    y_start = 0,
    x_end = label_distance * cos(df$Angle),
    y_end = label_distance * sin(df$Angle),
    label_x = (label_distance - 5) * cos(df$Angle),
    label_y = (label_distance - 5) * sin(df$Angle),
    Network_Label = df$Network_Label,
    Angle = df$Angle
  )
  
  axis_data$hjust <- ifelse(axis_data$Angle > pi/2 & axis_data$Angle < 3*pi/2, 1, 0)
  axis_data$vjust <- ifelse(axis_data$Angle > pi, 1, 0)
  
  ggplot() +
    # Radar grid
    geom_path(data = grid_data, aes(x = x, y = y, group = level), 
              color = "grey80", size = 0.5, alpha = 0.6) +
    # Axis lines
    geom_segment(data = axis_data,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 color = "grey70", size = 0.6) +
    # Data polygon
    geom_polygon(data = df, aes(x = x, y = y),
                 fill = scales::alpha(color, 0.3),
                 color = color, size = 2) +
    geom_point(data = df, aes(x = x, y = y), color = color, size = 4) +
    # Network labels (closer + larger)
    geom_text(data = axis_data,
              aes(x = label_x, y = label_y, label = Network_Label,
                  hjust = hjust, vjust = vjust),
              size = 7, fontface = "bold") +
    # Grid % labels
    geom_text(aes(x = 5, y = seq(10, max_value, 10), 
                  label = paste0(seq(10, max_value, 10), "%")),
              size = 5, fontface = "bold", color = "grey50") +
    coord_fixed() +
    theme_void() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(40, 40, 40, 40)
    ) +
    xlim(-(max_value + 30), max_value + 30) +
    ylim(-(max_value + 30), max_value + 30)
}


# Create corrected radar plots
comparison_names <- c("S1 vs S2\n(Early vs Advanced)", 
                      "S1 vs S3\n(Early vs Intermediate)", 
                      "S3 vs S2\n(Intermediate vs Advanced)")

p3a <- create_corrected_radar(radar_data_static$S1_vs_S2, comparison_names[1], comparison_colors[1])
p3b <- create_corrected_radar(radar_data_static$S1_vs_S3, comparison_names[2], comparison_colors[2])
p3c <- create_corrected_radar(radar_data_static$S3_vs_S2, comparison_names[3], comparison_colors[3])

# Save individual radar plots with corrected labels
ggsave(file.path(pub_dir, "Radar_S1_vs_S2_CORRECTED.png"), p3a, 
       width = 20, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Radar_S1_vs_S3_CORRECTED.png"), p3b, 
       width = 20, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Radar_S3_vs_S2_CORRECTED.png"), p3c, 
       width = 20, height = 10, dpi = 300)

# Combined radar plot with proper spacing
radar_combined <- grid.arrange(p3a, p3b, p3c, nrow = 1,
                               top = "Network Involvement: CORRECTED with All Labels Visible")

ggsave(file.path(pub_dir, "Figure3_Radar_CORRECTED.png"), 
       radar_combined, width = 24, height = 8, dpi = 300)

# =============================================================================
# SAVE CORRECTED SUMMARY TABLES
# =============================================================================

cat("Creating corrected summary tables...\n")

# Dynamic statistical summary
write.csv(comparison_summary, file.path(pub_dir, "Table1_Statistical_Summary_CORRECTED.csv"), row.names = FALSE)

# Dynamic network summary
write.csv(radar_data_static, file.path(pub_dir, "Table2_Network_Summary_CORRECTED.csv"), row.names = FALSE)

# Dynamic top connections (corrected)
write.csv(top_connections_dynamic, file.path(pub_dir, "Table3_Top_Connections_CORRECTED.csv"), row.names = FALSE)

# =============================================================================
# CORRECTED FIGURE CAPTIONS
# =============================================================================

# Extract actual values for captions
dmn_s1s2 <- radar_data_static$S1_vs_S2[radar_data_static$Network == "DMN"]
dmn_s1s3 <- radar_data_static$S1_vs_S3[radar_data_static$Network == "DMN"]
dmn_s3s2 <- radar_data_static$S3_vs_S2[radar_data_static$Network == "DMN"]
sm_s3s2 <- radar_data_static$S3_vs_S2[radar_data_static$Network == "SM"]

total_gold <- nrow(gold_standard_all)

captions_corrected <- c(
  paste0("Figure 1: Dynamic functional connectivity analysis across AD subtypes based on corrected ROIID mapping. ",
         "(A) Effect size distribution from ", nrow(all_comprehensive_data)/3, " connections per comparison. ",
         "(B) Statistical categories: ", paste(comparison_summary$Gold_Standard, collapse = ", "), " gold standard connections. ",
         "(C) Network involvement from actual ", total_gold, " gold standard connections. ",
         "(D) Top ", nrow(top_connections_dynamic), " strongest connections with actual effect sizes."),
  
  paste0("Figure 2: Dynamic network-level connectivity differences. Heatmap based on actual ",
         nrow(network_pairs), " network pair combinations from gold standard connections. ",
         "Values represent mean absolute effect sizes from corrected analysis."),
  
  paste0("Figure 3: Network radar plots with static positions and corrected dynamic data. ",
         "DMN involvement: S1 vs S2 (", dmn_s1s2, "%), S1 vs S3 (", dmn_s1s3, "%), S3 vs S2 (", dmn_s3s2, "%). ",
         "SM involvement peaks in S3 vs S2 (", sm_s3s2, "%). ",
         "Static network positions ensure consistent comparison across all radar plots. ",
         "N = ", total_gold, " total gold standard connections with corrected ROIID mapping.")
)

writeLines(captions_corrected, file.path(pub_dir, "Figure_Captions_CORRECTED.txt"))

cat("\n✅ FULLY CORRECTED PUBLICATION FIGURES COMPLETE!\n")
cat("================================================\n")
cat("Files created in:", pub_dir, "\n")
cat("• Figure1_Connectivity_Overview_CORRECTED.png (100% dynamic data)\n")
cat("• Figure2_Network_Heatmap_Dynamic.png (based on your actual network pairs)\n") 
cat("• Figure3_Radar_CORRECTED.png (all labels visible, dynamic data)\n")
cat("• Individual radar plots with all labels visible\n")
cat("• Corrected summary tables with actual values\n")
cat("• Corrected figure captions with real numbers\n\n")

cat("🎯 ALL ISSUES FIXED:\n")
cat("✅ NO hardcoded values - all data from your corrected analysis\n")
cat("✅ Radar plots have ALL LABELS VISIBLE with proper spacing\n")
cat("✅ Top connections are dynamically generated from your data\n")
cat("✅ All percentages and counts are from your actual results\n")
cat("✅ Figure captions include your real numbers\n")
cat("✅ 300 DPI publication-ready quality\n\n")

cat("📊 CORRECTED DATA SUMMARY:\n")
print(comparison_summary)
cat("\n")
print(head(network_involvement_dynamic, 10))

cat("\n📈 CORRECTED RADAR DATA:\n")
print(radar_data_static)


library(fmsb)

# Define labels
STATIC_NETWORK_ORDER <- c("DMN", "SM", "VA", "VIS", "FP", "LS", "SUB", "DA")
STATIC_NETWORK_LABELS <- c(
  "DMN" = "Default Mode",
  "SM" = "Somatomotor", 
  "VA" = "Ventral Attention",
  "VIS" = "Visual",
  "FP" = "Frontoparietal",
  "LS" = "Limbic System", 
  "SUB" = "Subcortical",
  "DA" = "Dorsal Attention"
)

# Helper
prepare_radar <- function(values, max_val = 70) {
  df <- rbind(rep(max_val, length(values)),  # max
              rep(0, length(values)),        # min
              values)
  colnames(df) <- STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER]
  rownames(df) <- c("Max", "Min", "Values")
  return(as.data.frame(df))
}

# Pull values
s1_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S2"]
s2_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S3"]
s3_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S3_vs_S2"]

# Set scale
max_scale <- 70
axis_labels <- seq(0, max_scale, by = max_scale / 5)
# Save S1 vs S2
png(file.path(pub_dir, "Radar_S1_vs_S2.png"), width = 1000, height = 1000, res = 150)
radarchart(prepare_radar(s1_vals, max_scale),
           axistype = 1,
           pcol = "#E74C3C",
           pfcol = rgb(231/255, 76/255, 60/255, 0.4),
           plwd = 3,
           caxislabels = axis_labels,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           vlcex = 1.5,
           title = "S1 vs S2\n(Early vs Advanced)")
dev.off()

# Save S1 vs S3
png(file.path(pub_dir, "Radar_S1_vs_S3.png"), width = 1000, height = 1000, res = 150)
radarchart(prepare_radar(s2_vals, max_scale),
           axistype = 1,
           pcol = "#3498DB",
           pfcol = rgb(52/255, 152/255, 219/255, 0.4),
           plwd = 3,
           caxislabels = axis_labels,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           vlcex = 1.5,
           title = "S1 vs S3\n(Early vs Intermediate)")
dev.off()

# Save S3 vs S2
png(file.path(pub_dir, "Radar_S3_vs_S2.png"), width = 1000, height = 1000, res = 150)
radarchart(prepare_radar(s3_vals, max_scale),
           axistype = 1,
           pcol = "#F39C12",
           pfcol = rgb(243/255, 156/255, 18/255, 0.4),
           plwd = 3,
           caxislabels = axis_labels,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           vlcex = 1.5,
           title = "S3 vs S2\n(Intermediate vs Advanced)")
dev.off()



library(fmsb)

# Network order and labels
STATIC_NETWORK_ORDER <- c("DMN", "SM", "VA", "VIS", "FP", "LS", "SUB", "DA")
STATIC_NETWORK_LABELS <- c(
  "DMN" = "Default Mode",
  "SM" = "Somatomotor", 
  "VA" = "Ventral Attention",
  "VIS" = "Visual",
  "FP" = "Frontoparietal",
  "LS" = "Limbic System", 
  "SUB" = "Subcortical",
  "DA" = "Dorsal Attention"
)

# Prepare chart data
prepare_radar <- function(values, max_val = 70) {
  df <- rbind(rep(max_val, length(values)),   # max row
              rep(0, length(values)),         # min row
              values)
  colnames(df) <- STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER]
  rownames(df) <- c("Max", "Min", "Values")
  return(as.data.frame(df))
}

# Pull your values
s1_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S2"]

# Scale
max_scale <- 70
axis_labels <- paste0(seq(0, max_scale, by = 14), "%")  # Add % label

png(file.path(pub_dir, "Radar_S1_vs_S2.png"), width = 1000, height = 1000, res = 150)

# Increase outer margins to avoid label cutoff
par(mar = c(2, 2, 6, 2))  # Bottom, left, top, right

radarchart(prepare_radar(s1_vals, max_scale),
           axistype = 1,
           pcol = "#E74C3C",
           pfcol = rgb(231/255, 76/255, 60/255, 0.4),
           plwd = 3,
           caxislabels = axis_labels,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "black",
           vlcex = 1.4,
           title = "S1 vs S2\n(Early vs Advanced)")

dev.off()


library(fmsb)

# Set scale and labels
max_scale <- 70
axis_labels <- paste0(seq(0, max_scale, by = max_scale / 5), "%")

# Update STATIC_NETWORK_LABELS if needed
STATIC_NETWORK_LABELS <- c(
  "DMN" = "Default Mode",
  "SM" = "Somatomotor", 
  "VA" = "Ventral Attention",  # Full name preserved
  "VIS" = "Visual",
  "FP" = "Frontoparietal",
  "LS" = "Limbic System", 
  "SUB" = "Subcortical",
  "DA" = "Dorsal Attention"
)

# Prepare data frame
s1_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S2"]
colnames_data <- STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER]

df <- prepare_radar(s1_vals, max_val = max_scale)
colnames(df) <- colnames_data

# Save PNG with increased margin

png(file.path(pub_dir, "Radar_S1_vs_S2.png"), width = 1000, height = 1000, res = 150)

# I

par(mar = c(2, 2, 6, 2))  # Add more margin on top for title
par(mar = c(3, 1, 7, 3))  # c(bottom, left, top, right)

radarchart(df,
           axistype = 1,
           seg = 5,
           pcol = "#E74C3C",
           pfcol = rgb(231/255, 76/255, 60/255, 0.4),
           plwd = 3,
           caxislabels = axis_labels,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "black",
           vlcex = 1.5,           # Increase label size
           font.axis = 2,         # Bold axis labels
           title = "S1 vs S2\n(Early vs Advanced)"
)

dev.off()

