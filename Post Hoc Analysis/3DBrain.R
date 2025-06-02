# Required libraries
library(ggseg3d)
library(dplyr)
library(plotly)
library(scales)
library(RColorBrewer)
library(tibble)
library(fs)
library(htmlwidgets)

# Access internal functions
get_atlas <- ggseg3d:::get_atlas
data_merge <- ggseg3d:::data_merge

# Set output directory
output_dir <- "C:/Users/theya/Downloads/ADNIFinal/FinalPlots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
pet_data <- read.csv("C:/Users/theya/Downloads/ADNIFinal/Brainplotk=3.csv")

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

# Function to process cluster data (same as original)
process_cluster_data <- function(pet_data, cluster_label) {
  cluster_raw <- pet_data %>%
    filter(PredictedLabel3 == cluster_label) %>%
    select(starts_with("V")) %>%
    summarise(across(everything(), mean)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("vertex") %>%
    rename(value = V1)
  
  cluster_data <- cluster_raw %>%
    left_join(vertex_to_atlas_mapping, by = "vertex") %>%
    arrange(desc(hemi))
  
  return(cluster_data)
}

# Process data for each cluster
cluster0_data <- process_cluster_data(pet_data, 0)
cluster1_data <- process_cluster_data(pet_data, 1)
cluster2_data <- process_cluster_data(pet_data, 2)

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
                   title = "Value",
                   ticklen = 2,
                   len = 0.5,
                   lenmode = "fraction",
                   y = 1,
                   yanchor = "top",
                   nticks = 10,
                   tickmode = "auto",
                   tickformat = ".2f",
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