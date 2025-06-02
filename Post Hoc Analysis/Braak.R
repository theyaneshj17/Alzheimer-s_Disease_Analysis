# Load libraries
library(dplyr)
library(readr)

# Load the data
df <- read_csv("C:/Users/theya/Downloads/ADNIFinal/FinalPlots3/Cluster_BrainRegion_ValuesFinal.csv")

library(stringr)


# Convert region names to lowercase for consistent matching
df <- df %>% mutate(region = tolower(region))

# Filter for Subtype S1 only
s1_data <- df %>% filter(cluster == "S3")

# Define region mapping for Braak stages
braak_map <- list(
  "Braak I–II" = c("entorhinal", "hippocampus"),
  "Braak III–IV" = c("amygdala", "parahippocampal", "fusiform", "lingual",
                     "insula", "inferior temporal", "middle temporal", "posterior cingulate", "inferior parietal"),
  "Braak V–VI" = c("lateral orbitofrontal", "superior temporal", "inferior frontal", "cuneus",
                   "caudal anterior cingulate", "rostral anterior cingulate", "pars triangularis",
                   "pars opercularis", "pars orbitalis", "caudal middle frontal", "rostral middle frontal",
                   "supramarginal", "lateral occipital", "precuneus", "superior parietal", "superior frontal",
                   "frontal pole", "paracentral", "postcentral", "precentral", "pericalcarine")
)

# ASCII bar plot function
ascii_bar <- function(value, max_val = NULL, bar_width = 20) {
  if (is.null(max_val)) max_val <- max(abs(value), na.rm = TRUE)
  bar_len <- round(abs(value) / max_val * bar_width)
  bar <- if (value >= 0) str_dup("▓", bar_len) else str_dup("░", bar_len)
  sprintf("%s %6.2f%%", bar, value * 100)
}

# Print ASCII bar chart for each Braak stage group
for (stage in names(braak_map)) {
  cat(sprintf("\n%s Tau Burden (Subtype S3)\n", stage))
  cat(strrep("=", nchar(stage) + 24), "\n")
  
  region_list <- braak_map[[stage]]
  
  plot_data <- s1_data %>%
    filter(region %in% region_list) %>%
    arrange(desc(avg_value))
  
  max_val <- max(abs(plot_data$avg_value), na.rm = TRUE)
  
  for (i in 1:nrow(plot_data)) {
    row <- plot_data[i, ]
    cat(sprintf("%-20s | %s\n", str_to_title(row$region), ascii_bar(row$avg_value, max_val)))
  }
  

}
