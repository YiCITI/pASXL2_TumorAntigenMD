# ==============================================================================
# Script Name: RMSF_Analysis.R
# Description: 
#   1. Reads RMSF (Root Mean Square Fluctuation) data from a CSV file.
#   2. Reshapes data from Wide format (columns: 1_pos, 1_neg...) to Long format.
#   3. Visualizes RMSF distribution using grouped Boxplots (Positive vs Negative).
#   4. Performs Wilcoxon tests between groups for each residue.
#
# Input Data Requirement:
#   Format: .csv (Header = TRUE)
#   Column Names Pattern: "{ResidueID}_{Group}" (e.g., "1_pos", "1_neg", "2_pos")
#   Path: ./data/data.csv
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(data.table) 

# ==============================================================================
# 1. Configuration & Data Loading
# ==============================================================================

# [Config]: Paths (Relative paths for GitHub portability)
input_file  <- "./data/data.csv"
output_dir  <- "./results"
output_file <- file.path(output_dir, "RMSF_Boxplot.pdf")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# [Config]: Visualization Colors
# pos (e.g., Tumor/Mutant) -> Red, neg (e.g., Normal/WT) -> Blue
cols <- c("pos" = "#E31A1C", "neg" = "#377EB8")

message(paste("Reading data from:", input_file))

# Check if file exists
if(!file.exists(input_file)){
  stop(paste("File not found:", input_file, "\nPlease ensure your data is placed in the 'data' folder."))
}

# Read data (auto-detect separator)
df_raw <- fread(input_file, header = TRUE, check.names = FALSE)

# ==============================================================================
# 2. Data Processing (Wide -> Long)
# ==============================================================================

message("Processing data structure...")

plot_df <- df_raw %>%
  # 1. Pivot to Long format: Convert all columns to Key-Value pairs
  pivot_longer(cols = everything(), names_to = "RawName", values_to = "Value") %>%
  
  # 2. Split RawName into ID and Group based on the underscore '_'
  mutate(
    # Extract part before "_" as Residue ID
    ID = sub("_.*", "", RawName),
    # Extract part after "_" as Group (pos/neg)
    Group = sub(".*_", "", RawName)
  )

# 3. Set Factor Levels
# Sort ID numerically (1, 2, ... 9, 10) instead of alphabetically (1, 10, 2...)
unique_ids <- sort(as.numeric(unique(plot_df$ID)))
plot_df$ID <- factor(plot_df$ID, levels = as.character(unique_ids))

# Set Group order: neg (Blue) first, pos (Red) second
plot_df$Group <- factor(plot_df$Group, levels = c("neg", "pos"))

# Remove any NA values introduced during conversion
plot_df <- plot_df %>% filter(!is.na(Value))

message(paste("Data ready. Total observations:", nrow(plot_df)))

# ==============================================================================
# 3. Visualization
# ==============================================================================

p <- ggplot(plot_df, aes(x = ID, y = Value, fill = Group, color = Group)) +
  
  # 1. Boxplot Body (fatten = 0 hides the default black median line)
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,
    fatten = 0,                
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.shape = NA
  ) +
  
  # 2. Custom White Median Line
  stat_summary(
    fun = median,
    geom = "errorbar",        
    aes(ymax = ..y.., ymin = ..y.., group = Group), 
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,             
    color = "white",         
    size = 0.8               
  ) +
  
  # 3. Colors
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  # 4. Theme
  theme_classic() +
  
  # 5. Statistical Comparison (Wilcoxon test within each Residue ID)
  stat_compare_means(
    aes(group = Group), 
    label = "p.signif", 
    method = "wilcox.test",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    ),
    hide.ns = FALSE,
    vjust = -0.5,
    color = "black"
  ) + 
  
  # 6. Axis adjustments (Expand top margin for stars)
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  
  labs(x = "Residue ID", y = "RMSF (nm)", title = "") +
  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black"),
    panel.border = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ==============================================================================
# 4. Save Output
# ==============================================================================
ggsave(output_file, p, width = 8, height = 6)
message(paste("Plot saved to:", output_file))



































