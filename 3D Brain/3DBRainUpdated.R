# Required libraries
library(ggseg3d)
library(dplyr)
library(plotly)
library(scales)
library(RColorBrewer)
library(tibble)
library(fs)
library(htmlwidgets)

pet_data <- read.csv("C:/Users/theya/Downloads/ADNIFinal/Brainplotk=3.csv")
library(tidyr)
# Access internal functions
get_atlas <- ggseg3d:::get_atlas
data_merge <- ggseg3d:::data_merge


output_dir <- "C:/Users/theya/Downloads/ADNIFinal/FinalPlots3"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

vertex_cols <- paste0("V", 1:68)

global_min <- min(pet_data[, vertex_cols], na.rm = TRUE)

shift_amount <- abs(global_min) + 0.001

vertex_data <- pet_data[, vertex_cols] + shift_amount


pet_data$total_tau <- rowSums(vertex_data, na.rm = TRUE)

normalized_data <- pet_data

normalized_data[, vertex_cols] <- sweep(vertex_data, 1, pet_data$total_tau, FUN = "/")


normalized_data$SubjectID <- pet_data$SubjectID

normalized_data$PredictedLabel3 <- pet_data$PredictedLabel3

# Compute average regional proportion within each predicted subtype group (S1/S2/S3)
library(dplyr)

group_avg <- normalized_data %>%
  group_by(PredictedLabel3) %>%
  summarise(across(all_of(vertex_cols), mean, na.rm = TRUE))

# Optional: Save group-averaged normalized data
write.csv(group_avg, file.path(output_dir, "GroupAvg_NormalizedTau.csv"), row.names = FALSE)

# View group-average normalized values
print(group_avg)
check_sums <- group_avg %>%
  arrange(PredictedLabel3) %>%
  rowwise() %>%
  mutate(sum_props = sum(c_across(all_of(vertex_cols)))) %>%
  select(PredictedLabel3, sum_props)

print(check_sums)


print(round(check_sums, 4))  # Should print 1.0000 or very close
table(normalized_data$PredictedLabel3)
summary(pet_data$total_tau)


# Create vertex to atlas mapping (keeping the same mapping as original)
vertex_to_atlas_mapping <- data.frame(
  vertex = paste0("V", 1:68),
  region_name = c(
    # Right hemi (V1-V34)
    "lateral orbitofrontal", "pars orbitalis", "frontal pole",
    "medial orbitofrontal", "pars triangularis", "pars opercularis",
    "rostral middle frontal", "superior frontal", "caudal middle frontal",
    "precentral", "paracentral", "rostral anterior cingulate",
    "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate",
    "postcentral", "supramarginal", "superior parietal",
    "inferior parietal", "precuneus", "cuneus", "pericalcarine",
    "lateral occipital", "lingual", "fusiform", "parahippocampal",
    "entorhinal", "temporal pole", "inferior temporal", "middle temporal",
    "bankssts", "superior temporal", "transverse temporal", "insula",
    # Left hemi (V35-V68, same order)
    "lateral orbitofrontal", "pars orbitalis", "frontal pole",
    "medial orbitofrontal", "pars triangularis", "pars opercularis",
    "rostral middle frontal", "superior frontal", "caudal middle frontal",
    "precentral", "paracentral", "rostral anterior cingulate",
    "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate",
    "postcentral", "supramarginal", "superior parietal",
    "inferior parietal", "precuneus", "cuneus", "pericalcarine",
    "lateral occipital", "lingual", "fusiform", "parahippocampal",
    "entorhinal", "temporal pole", "inferior temporal", "middle temporal",
    "bankssts", "superior temporal", "transverse temporal", "insula"
  ),
  hemi = rep(c("right", "left"), each = 34)
)

# Define specific views we want
views <- list(
  lateral_right = list(x = 2, y = 0, z = 0),    # Right side view
  lateral_left = list(x = -2, y = 0, z = 0),    # Left side view
  superior = list(x = 0, y = 0, z = 2)          # Top view
)


process_cluster_data <- function(group_avg_row, cluster_label) {
  # Extract 1 row (cluster-specific averaged values)
  df <- group_avg_row[group_avg_row$PredictedLabel3 == cluster_label, ]
  
  # Convert wide to long
  long_df <- df %>%
    select(all_of(vertex_cols)) %>%
    pivot_longer(cols = everything(), names_to = "vertex", values_to = "value")
  
  # Join with atlas
  long_df <- long_df %>%
    left_join(vertex_to_atlas_mapping, by = "vertex") %>%
    arrange(desc(hemi))
  
  return(long_df)
}




# Process data for each cluster
cluster0_data <- process_cluster_data(group_avg, 0)
cluster1_data <- process_cluster_data(group_avg, 1)
cluster2_data <- process_cluster_data(group_avg, 2)

# Function to summarize brain data
summarize_brain_data <- function(cluster_data) {
  cluster_data %>%
    group_by(region_name, hemi) %>%
    summarise(
      value = mean(value),
      .groups = "drop"
    ) %>%
    rename(region = region_name)
}

# Create brain data for each cluster
brain_data_0 <- summarize_brain_data(cluster0_data)
brain_data_1 <- summarize_brain_data(cluster1_data)
brain_data_2 <- summarize_brain_data(cluster2_data)


# Add cluster label to each data frame
brain_data_0$cluster <- "S1"
brain_data_1$cluster <- "S2"
brain_data_2$cluster <- "S3"


brain_data_avg_0 <- brain_data_0 %>%
  group_by(region) %>%
  summarise(avg_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(cluster = "S1") %>%
  arrange(desc(avg_value))

brain_data_avg_1 <- brain_data_1 %>%
  group_by(region) %>%
  summarise(avg_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(cluster = "S2") %>%
  arrange(desc(avg_value))

brain_data_avg_2 <- brain_data_2 %>%
  group_by(region) %>%
  summarise(avg_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(cluster = "S3") %>%
  arrange(desc(avg_value))

# Combine all cluster data into one dataframe
combined_brain_data <- bind_rows(brain_data_avg_0, brain_data_avg_1, brain_data_avg_2)

# Save to CSV
write.csv(combined_brain_data,
          file = file.path(output_dir, "Cluster_BrainRegion_ValuesFinal.csv"),
          row.names = FALSE)



# Modified function to create brain visualization
create_brain_viz <- function(brain_data, view_name, view) {
  # Calculate value range across all clusters
  value_range <- c(
    min(c(brain_data_0$value, brain_data_1$value, brain_data_2$value)),
    max(c(brain_data_0$value, brain_data_1$value, brain_data_2$value))
  )
  
  # Create color palette
  colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
  
  # Create color mapping function
  color_func <- gradient_n_pal(
    colours = colors,
    values = seq(value_range[1], value_range[2], length.out = length(colors)),
    space = "Lab"
  )
  
  # Get atlas data and merge
  atlas3d <- get_atlas(dk_3d, surface = "inflated", hemisphere = c("left", "right"))
  atlas3d_unique <- atlas3d[!duplicated(paste(atlas3d$hemi, atlas3d$region)),]
  atlas3d <- data_merge(brain_data, atlas3d_unique)
  
  # Apply color mapping
  atlas3d$new_col <- color_func(atlas3d$value)
  
  # Create plot
  p <- plot_ly() %>%
    layout(
      title = paste("Brain Visualization -", view_name),
      scene = list(
        xaxis = list(visible = FALSE, showgrid = FALSE),
        yaxis = list(visible = FALSE, showgrid = FALSE),
        zaxis = list(visible = FALSE, showgrid = FALSE),
        camera = list(eye = view),
        bgcolor = "white"
      )
    )
  
  # Add brain mesh traces
  for (tt in 1:nrow(atlas3d)) {
    col <- rep(atlas3d$new_col[tt], nrow(atlas3d$mesh[[tt]]$faces))
    col <- ifelse(is.na(col), "darkgrey", col)
    
    p <- add_trace(p, 
                   x = atlas3d$mesh[[tt]]$vertices$x,
                   y = atlas3d$mesh[[tt]]$vertices$y,
                   z = atlas3d$mesh[[tt]]$vertices$z,
                   i = atlas3d$mesh[[tt]]$faces$i - 1,
                   j = atlas3d$mesh[[tt]]$faces$j - 1,
                   k = atlas3d$mesh[[tt]]$faces$k - 1,
                   facecolor = col,
                   type = "mesh3d",
                   showscale = FALSE,
                   name = atlas3d$region[tt]
    )
  }
  
  # Add colorbar
  p <- add_trace(p,
                 x = 0, y = 0, z = 0,
                 intensity = seq(value_range[1], value_range[2], length.out = length(colors)),
                 colorscale = list(seq(0, 1, length.out = length(colors)), colors),
                 type = "mesh3d",
                 colorbar = list(
                   title = list(
                     text = "Value",
                     font = list(
                       size = 35,
                       family = "Arial Black"  # Bold title
                     )
                   ),
                   tickfont = list(
                     size = 25,
                     family = "Arial Black"  # Bold tick labels
                   ),
                   ticklen = 2,
                   len = 0.5,
                   lenmode = "fraction",
                   y = 1,
                   yanchor = "top",
                   nticks = 10,
                   tickmode = "auto",
                   tickformat = ".3f",
                   thickness = 20,
                   outlinewidth = 1
                 )
  )
  
  
  return(p)
}

# Function to save views for a cluster
save_cluster_views <- function(brain_data, cluster_num) {
  # Create cluster-specific directory
  cluster_dir <- file.path(output_dir, paste0("cluster", cluster_num))
  dir.create(cluster_dir, showWarnings = FALSE)
  
  # Save each view
  for (view_name in names(views)) {
    view_viz <- create_brain_viz(brain_data, view_name, views[[view_name]])
    htmlwidgets::saveWidget(
      view_viz,
      file.path(cluster_dir, paste0("cluster", cluster_num, "_", view_name, ".html")),
      selfcontained = TRUE
    )
  }
}

# Save visualizations for all clusters
save_cluster_views(brain_data_0, 0)
save_cluster_views(brain_data_1, 1)
save_cluster_views(brain_data_2, 2)

