library(fmsb)

# Define static network order and labels
STATIC_NETWORK_ORDER <- c("DMN", "SM", "VA", "VIS", "FP", "LS", "SUB", "DA")
STATIC_NETWORK_LABELS <- c(
  "DMN" = "DMN", "SM" = "SM", "VA" = "VA", "VIS" = "VIS",
  "FP" = "FP", "LS" = "LS", "SUB" = "SUB", "DA" = "DA"
)

# Static configuration parameters
MAX_SCALE <- 70
NUM_RINGS <- 5
AXIS_LABELS <- c("0%", "14%", "28%", "42%", "56%", "70%")

# Helper function to prepare radar data
prepare_radar <- function(values, max_val = MAX_SCALE) {
  if (length(values) != length(STATIC_NETWORK_ORDER)) {
    stop("Values must match the number of networks")
  }
  
  # Create the radar data frame structure
  df <- rbind(rep(max_val, length(values)),  # max values
              rep(0, length(values)),        # min values  
              values)                        # actual values
  
  # Use simple network labels
  network_labels <- STATIC_NETWORK_LABELS[STATIC_NETWORK_ORDER]
  
  # Set proper column names
  colnames(df) <- network_labels
  rownames(df) <- c("Max", "Min", "Values")
  
  return(as.data.frame(df))
}

# Extract values in the correct network order
s1_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S2"]
s2_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S1_vs_S3"]
s3_vals <- radar_data_static[STATIC_NETWORK_ORDER, "S3_vs_S2"]

# Create combined 3-panel radar chart with enhanced formatting
create_combined_radar_chart <- function() {
  # Create large PNG for 3 panels side by side
  png(file.path(pub_dir, "Combined_Radar_Charts.png"), 
      width = 6000, height = 2200, res = 200)
  
  # Set up 3 panels in 1 row with very tight margins to bring labels closer
  par(mfrow = c(1, 3), mar = c(0.2, 0.2, 1.5, 0.2), oma = c(0, 0, 0, 0))
  
  # Set text formatting to bold
  par(font.axis = 2, font.lab = 2, font.main = 2)  # Make all text bold
  
  # Panel 1: S1 vs S2
  radarchart(prepare_radar(s1_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.2,
             pcol = "#E74C3C",
             pfcol = rgb(231/255, 76/255, 60/255, 0.4),
             plwd = 4,
             axislabcol = "black",
             vlcex = 4.5,          # Even larger labels
             calcex = 5.0,         # Even larger percentages
             centerzero = TRUE,
             maxmin = FALSE)       # This can help with label positioning
  
  # Panel 2: S1 vs S3
  radarchart(prepare_radar(s2_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.2,
             pcol = "#3498DB",
             pfcol = rgb(52/255, 152/255, 219/255, 0.4),
             plwd = 4,
             axislabcol = "black",
             vlcex = 4.5,          # Even larger labels
             calcex = 5.0,         # Even larger percentages
             centerzero = TRUE,
             maxmin = FALSE)       # This can help with label positioning
  
  # Panel 3: S3 vs S2
  radarchart(prepare_radar(s3_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.2,
             pcol = "#F39C12",
             pfcol = rgb(243/255, 156/255, 18/255, 0.4),
             plwd = 4,
             axislabcol = "black",
             vlcex = 4.5,          # Even larger labels
             calcex = 5.0,         # Even larger percentages
             centerzero = TRUE,
             maxmin = FALSE)       # This can help with label positioning
  
  dev.off()
}

# Alternative approach with custom text rendering for maximum boldness
create_enhanced_radar_chart <- function() {
  # Create even larger PNG for maximum quality
  png(file.path(pub_dir, "Enhanced_Radar_Charts.png"), 
      width = 8000, height = 3000, res = 250)
  
  # Set up 3 panels with minimal margins
  par(mfrow = c(1, 3), mar = c(0.1, 0.1, 1, 0.1), oma = c(0, 0, 0, 0))
  
  # Force all text to be bold and larger
  par(font = 2, font.axis = 2, font.lab = 2, font.main = 2, cex = 1.2)
  
  # Panel 1: S1 vs S2
  radarchart(prepare_radar(s1_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.5,
             pcol = "#E74C3C",
             pfcol = rgb(231/255, 76/255, 60/255, 0.4),
             plwd = 5,
             axislabcol = "black",
             vlcex = 5.5,          # Maximum label size
             calcex = 5.5,         # Maximum percentage size
             centerzero = TRUE,
             maxmin = FALSE)
  
  # Panel 2: S1 vs S3
  radarchart(prepare_radar(s2_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.5,
             pcol = "#3498DB",
             pfcol = rgb(52/255, 152/255, 219/255, 0.4),
             plwd = 5,
             axislabcol = "black",
             vlcex = 5.5,          # Maximum label size
             calcex = 5.5,         # Maximum percentage size
             centerzero = TRUE,
             maxmin = FALSE)
  
  # Panel 3: S3 vs S2
  radarchart(prepare_radar(s3_vals, MAX_SCALE),
             axistype = 1,
             seg = NUM_RINGS,
             caxislabels = AXIS_LABELS,
             cglcol = "grey60",
             cglty = 1,
             cglwd = 1.5,
             pcol = "#F39C12",
             pfcol = rgb(243/255, 156/255, 18/255, 0.4),
             plwd = 5,
             axislabcol = "black",
             vlcex = 5.5,          # Maximum label size
             calcex = 5.5,         # Maximum percentage size
             centerzero = TRUE,
             maxmin = FALSE)
  
  dev.off()
}

# Run both versions
create_combined_radar_chart()    # Standard version with bold formatting
create_enhanced_radar_chart()    # Enhanced version with maximum boldness

# Print configuration summary
cat("Enhanced Radar Chart Configuration:\n")
cat("- Standard version: 6000x2200 pixels with bold formatting\n")
cat("- Enhanced version: 8000x3000 pixels with maximum boldness\n")
cat("- Ultra-tight margins: mar = c(0.1, 0.1, 1, 0.1) for closest labels\n")
cat("- Bold text forced with: font = 2, font.axis = 2, font.lab = 2\n")
cat("- Maximum label sizes: vlcex = 5.5, calcex = 5.5\n")
cat("- Networks (order):", paste(STATIC_NETWORK_ORDER, collapse = ", "), "\n")  
cat("- Scale: 0 to", MAX_SCALE, "with bold percentage labels\n")
cat("- Number of rings:", NUM_RINGS, "\n") 
cat("- Axis labels:", paste(AXIS_LABELS, collapse = ", "), "\n")