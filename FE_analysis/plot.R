# ==============================================================================
# Script Name: Binding_Energy_Analysis.R
# Description: 
#   1. Reads binding energy data (CSV) for different simulation groups.
#   2. Merges data and performs data cleaning.
#   3. Visualizes Binding Free Energy distributions using Boxplots.
#   4. Performs statistical comparisons (Wilcoxon test) with dynamic positioning.
#
# Input Data Requirement:
#   Format: .csv (Header = TRUE)
#   Column: Must contain a column named "DELTA TOTAL" (Binding Energy).
#   Path: ./data/{AnalysisName}/Phos.csv, Non_phos.csv, etc.
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)

# ==============================================================================
# 1. Configuration [Modify this section]
# ==============================================================================

# [Config]: Analysis Name (Folder name in ./data and Output filename prefix)
# Change this to "Figure3A" or "Figure3B" as needed
analysis_name <- "Figure3B"

# [Config]: Directories (Relative paths for GitHub portability)
base_data_dir <- "./data"
output_dir <- "./results"

# Construct specific paths
work_dir <- file.path(base_data_dir, analysis_name)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# [Config]: File Names and corresponding Group Labels
# Ensure these lists have the same length and order
file_names <- c("Phos.csv", "Non_phos.csv", "Prot_phos.csv", "Muted_phos1.csv", "Muted_phos2.csv")
group_names <- c("Phos", "Non_phos", "Prot_phos", "Muted_phos1", "Muted_phos2")

# [Config]: Visualization Colors
my_colors <- c(
  "Phos"        = "#C0392B", 
  "Non_phos"    = "#2874A6", 
  "Prot_phos"   = "#117A65", 
  "Muted_phos1" = "#D35400", 
  "Muted_phos2" = "#B7950B"  
)

# [Config]: Statistical Comparisons
# Define pairs for Wilcoxon test
my_comparisons <- list(
  c("Phos", "Non_phos"),      
  c("Non_phos", "Prot_phos"), 
  c("Phos", "Prot_phos")      
)

# ==============================================================================
# 2. Data Loading & Processing
# ==============================================================================
message(paste("Loading data from:", work_dir))

data_list <- list()

for (i in seq_along(file_names)) {
  file_path <- file.path(work_dir, file_names[i])
  
  if(file.exists(file_path)){
    temp_df <- fread(file_path, header = TRUE)
    
    # Ensure necessary column exists
    if (!"DELTA TOTAL" %in% colnames(temp_df)) {
      warning(paste("Column 'DELTA TOTAL' not found in", file_names[i]))
      next
    }
    
    temp_df$Group <- group_names[i]
    data_list[[i]] <- temp_df
  } else {
    warning(paste("File not found:", file_path))
  }
}

if (length(data_list) == 0) stop("No valid data loaded. Check your data directory.")

# Merge Data
plot_data <- do.call(rbind, data_list)

# Data Cleaning
# Convert Group to Factor with specified levels order
plot_data$Group <- factor(plot_data$Group, levels = group_names)
plot_data$`DELTA TOTAL` <- as.numeric(plot_data$`DELTA TOTAL`)

# Remove NAs
plot_data <- plot_data %>% filter(!is.na(`DELTA TOTAL`))

message(paste("Total samples loaded:", nrow(plot_data)))

# ==============================================================================
# 3. Dynamic Calculation for Plot Aesthetics
# ==============================================================================

# Calculate Y-axis limits to determine P-value bracket positions
max_y <- max(plot_data$`DELTA TOTAL`, na.rm = TRUE)
min_y <- min(plot_data$`DELTA TOTAL`, na.rm = TRUE)
y_range <- max_y - min_y
step_size <- y_range * 0.05 

# Define Y positions for the comparison brackets manually but based on data max
# 1. First level comparisons
# 2. Second level comparisons (offset to avoid overlap)
# 3. Top level comparisons
y_positions <- c(
  max_y + step_size,         # Phos vs Non_phos
  max_y + step_size * 2.5,   # Non_phos vs Prot_phos
  max_y + step_size * 5.0    # Phos vs Prot_phos
)

# ==============================================================================
# 4. Plotting
# ==============================================================================

p <- ggplot(plot_data, aes(x = Group, y = `DELTA TOTAL`, fill = Group)) +
  
  # Boxplot layer
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) + 
  
  # Custom Colors
  scale_fill_manual(values = my_colors) +
  
  # Statistical Comparisons (Wilcoxon Test)
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    label.y = y_positions,   # Apply dynamic Y positions
    tip.length = 0.01,       # Short bracket tips for cleaner look
    size = 5,                # Star size
    hide.ns = FALSE,
    color = "black",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  
  # Theme and Labels
  theme_classic() + 
  
  labs(x = "", 
       y = "Binding Free Energy (kcal/mol)", 
       title = paste(analysis_name, "- Binding Energy Comparison")) +
  
  # Custom Theme Adjustments
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 0), 
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "none" 
  )

# Add top margin to accommodate high p-value brackets
p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

# ==============================================================================
# 5. Saving Results
# ==============================================================================
output_filename <- file.path(output_dir, paste0(analysis_name, "_Binding_Energy_Boxplot.pdf"))
ggsave(output_filename, p, width = 7, height = 6)
message(paste("Plot saved to:", output_filename))


