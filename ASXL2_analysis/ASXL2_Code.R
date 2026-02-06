# ==============================================================================
# Script Name: Pancancer_Protein_Abundance_Analysis.R
# Description: 
#   1. Extracts specific gene protein abundance from CPTAC data across multiple cancers.
#   2. Calculates Signed Meta P-value using Stouffer's Z-score method.
#   3. Generates a boxplot with statistical annotations.
#
# Input Data Structure Requirement:
#   ./data/{CancerType}/{CancerType}_proteomics_gene_abundance_log2_reference_intensity_normalized_{Tumor/Normal}.txt
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)      
library(org.Hs.eg.db) 

# ==============================================================================
# 1. Configuration & Preparation
# ==============================================================================

# Target gene symbol
target_gene <- "ASXL2"

# [Config]: Define base paths (Relative paths recommended for GitHub)
# Ensure your local data structure matches this path
base_dir <- "./data/CPTACPancancer" 
output_dir <- "./results"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# List of cancer types to analyze
cancer_types <- c("LUAD", "LSCC", "PDAC", "HNSCC", "CCRCC", "UCEC",    
                  "OV", "COAD", "GBM", "BRCA")

# 1.1 Retrieve Ensembl ID
message(paste("Fetching Ensembl ID for:", target_gene))
target_ensembl <- mapIds(org.Hs.eg.db, 
                         keys = target_gene, 
                         column = "ENSEMBL", 
                         keytype = "SYMBOL", 
                         multiVals = "first")

if (is.na(target_ensembl)) stop(paste("Ensembl ID not found for", target_gene))
message(paste("Target Gene ID:", target_ensembl))

# 1.3 Initialize list to store data
plot_data_list <- list()

# ==============================================================================
# 2. Helper Function: Extract Data from Single File
# ==============================================================================
extract_data <- function(filepath, cancer_name, group_type) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  dt <- fread(filepath, header = TRUE, sep = "\t")
  
  # Remove version number from Ensembl ID (e.g., ENSG0000.1 -> ENSG0000)
  dt$clean_id <- sub("\\.[0-9]+$", "", dt[[1]])
  target_rows <- dt[clean_id == target_ensembl]
  
  if (nrow(target_rows) == 0) return(NULL)
  
  # Extract numeric columns (samples)
  numeric_cols <- colnames(dt)[2:ncol(dt)] 
  mat <- as.matrix(target_rows[, ..numeric_cols])
  
  # Handle cases with multiple probes/rows (take median) or single row
  if (nrow(mat) > 1) {
    sample_values <- apply(mat, 2, median, na.rm = TRUE)
  } else {
    sample_values <- as.numeric(mat[1, ])
    names(sample_values) <- colnames(mat)
  }
  
  df_res <- data.frame(
    Sample = names(sample_values),
    Value = sample_values,
    Group = group_type,
    Cancer = cancer_name,
    stringsAsFactors = FALSE
  )
  return(df_res)
}

# ==============================================================================
# 3. Process All Cancer Types
# ==============================================================================

for (cancer in cancer_types) {
  message(paste("Processing:", cancer, "..."))
  dir_path <- file.path(base_dir, cancer)
  
  # Construct file paths (Ensure filename format matches your local data)
  f_tumor <- file.path(dir_path, paste0(cancer, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))
  f_normal <- file.path(dir_path, paste0(cancer, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"))
  
  dat_t <- extract_data(f_tumor, cancer, "Tumor")
  dat_n <- extract_data(f_normal, cancer, "Normal")
  
  if (!is.null(dat_t)) plot_data_list[[paste0(cancer, "_T")]] <- dat_t
  if (!is.null(dat_n)) plot_data_list[[paste0(cancer, "_N")]] <- dat_n
}

# ==============================================================================
# 4. Data Merging & Validation
# ==============================================================================

if (length(plot_data_list) == 0) stop("No valid data found. Please check paths and filenames.")

final_df <- bind_rows(plot_data_list)
final_df <- final_df %>% filter(!is.na(Value))
final_df$Group <- factor(final_df$Group, levels = c("Normal", "Tumor"))

# Update cancer levels based on available data
available_cancers <- unique(final_df$Cancer)
final_levels <- intersect(cancer_types, available_cancers)
final_df$Cancer <- factor(final_df$Cancer, levels = final_levels)

# ==============================================================================
# 4.5 Calculate Signed Meta P-value (Stouffer's Z-score Method)
# ==============================================================================
message("Calculating Signed Meta P-value (Stouffer's Method)...")

z_scores <- c()

for (c in final_levels) {
  sub_df <- final_df %>% filter(Cancer == c)
  
  # Calculation requires both Tumor and Normal groups
  if (n_distinct(sub_df$Group) == 2) {
    
    # 1. Wilcoxon test (two-sided)
    res <- wilcox.test(Value ~ Group, data = sub_df)
    p_val <- res$p.value
    
    # Handle extremely small p-values to avoid Inf Z-scores
    if (p_val < 1e-300) p_val <- 1e-300
    
    # 2. Determine direction (Tumor median vs Normal median)
    med_t <- median(sub_df$Value[sub_df$Group == "Tumor"], na.rm=TRUE)
    med_n <- median(sub_df$Value[sub_df$Group == "Normal"], na.rm=TRUE)
    
    # Sign: +1 if Tumor > Normal, -1 if Tumor < Normal
    direction <- sign(med_t - med_n)
    if(direction == 0) direction <- 1 
    
    # 3. Convert to signed Z-score (Inverse Normal Transform)
    z <- qnorm(1 - p_val / 2) * direction
    
    z_scores <- c(z_scores, z)
  }
}

# 4. Combine Z-scores using Stouffer's formula
if (length(z_scores) > 0) {
  k <- length(z_scores)
  meta_z <- sum(z_scores) / sqrt(k)
  
  # 5. Convert Meta Z back to P-value (two-sided)
  meta_p_val <- 2 * (1 - pnorm(abs(meta_z)))
  
  # 6. Generate label text
  sign_symbol <- ifelse(meta_z > 0, "(+)", "(-)")
  
  if (meta_p_val < 0.0001) {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " < 0.0001")
  } else {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " = ", sprintf("%.4f", meta_p_val))
  }
} else {
  meta_label <- "Meta P: N/A"
}

message(paste("Calculation Complete:", meta_label))

# ==============================================================================
# 5. Visualization
# ==============================================================================

cols <- c("Tumor" = "#E31A1C", "Normal" = "#377EB8")

p <- ggplot(final_df, aes(x = Cancer, y = Value, fill = Group, color = Group)) +
  
  # Boxplot layer
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,
    fatten = 0,             
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.colour = NULL,
    outlier.shape = NA
  ) +
  
  # Median line layer (White)
  stat_summary(
    fun = median,
    geom = "errorbar",      
    aes(ymax = ..y.., ymin = ..y.., group = Group), 
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,            
    color = "white",        
    size = 0.8              
  ) +
  
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  theme_classic() +
  
  # Statistical comparison (Wilcoxon test)
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", 
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     hide.ns = FALSE,
                     vjust = -0.5,
                     color = "black") + 
  
  # Add Meta P-value annotation
  annotate("text", 
           x = 0.6,                          
           y = max(final_df$Value) * 1.05,  
           label = meta_label, 
           hjust = 0,                        
           size = 5,                         
           fontface = "bold",                
           color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  labs(x = "", y = "Log2 Intensity", 
       title = paste0(target_gene, " Abundance (Proteomics)")) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
    panel.border = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save plot
output_filename <- file.path(output_dir, paste0(target_gene, "_Pancancer_Protein_Boxplot.pdf"))
ggsave(output_filename, p, width = 7, height = 6)

message(paste("Plot saved to:", output_filename))


























# ==============================================================================
# Script Name: Pancancer_mRNA_Expression_Analysis.R
# Description: 
#   1. Extracts specific gene mRNA expression (TPM) from CPTAC Excel files.
#   2. Maps Gene Symbols to Ensembl IDs.
#   3. Calculates Signed Meta P-value using Stouffer's Z-score method.
#   4. Generates a boxplot with statistical annotations (White median line style).
#
# Input Data Requirement:
#   Format: .xlsx (Excel)
#   Path: ./data/{CancerType}/{CancerType}_TPM.xlsx
#   Structure: First column contains Ensembl IDs (e.g., ENSG0000xxx.y).
#              Sample names ending in 'A' are Normal, 'T' are Tumor.
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)       # For statistical tests
library(org.Hs.eg.db) # For ID conversion
library(readxl)       # For reading .xlsx files

# ==============================================================================
# 1. Configuration & Preparation
# ==============================================================================

target_gene <- "ASXL2"

# [Config]: Define base paths (Relative paths recommended for GitHub)
base_dir <- "./data/CPTACPancancer"
output_dir <- "./results"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1.1 Retrieve Ensembl ID
message(paste("Fetching Ensembl ID for:", target_gene))
target_ensembl <- mapIds(org.Hs.eg.db, 
                         keys = target_gene, 
                         column = "ENSEMBL", 
                         keytype = "SYMBOL", 
                         multiVals = "first")

if (is.na(target_ensembl)) stop(paste("Ensembl ID not found for", target_gene))
message(paste("Target Gene ID:", target_ensembl))

# 1.2 Define Cancer Types
cancer_types <- c("LUAD", "LSCC", "PDAC", "HNSCC", "CCRCC", "UCEC",    
                  "OV", "COAD", "GBM", "BRCA")

# Define cancers that have paired Normal samples (Mixed in the same file)
cancers_with_normal <- c("CCRCC", "HNSCC", "LSCC", "LUAD", "PDAC", "UCEC")

# 1.3 Initialize list to store data
plot_data_list <- list()

# ==============================================================================
# 2. Helper Function: Extract RNA Data (Excel Format)
# ==============================================================================
extract_rna_data <- function(cancer_name) {
  
  # Construct file path
  file_path <- file.path(base_dir, cancer_name, paste0(cancer_name, "_TPM.xlsx"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  message(paste("Reading:", cancer_name, "..."))
  
  # Read Excel file (suppress messages to keep console clean)
  df <- suppressMessages(read_excel(file_path))
  df <- as.data.frame(df)
  
  # Clean IDs (Remove version number, e.g., ENSG0000.1 -> ENSG0000)
  raw_ids <- df[[1]]
  clean_ids <- sub("\\.[0-9]+$", "", raw_ids)
  
  # Find target gene index
  target_idx <- which(clean_ids == target_ensembl)
  
  if (length(target_idx) == 0) {
    message(paste("Skip:", cancer_name, "- Target gene not found"))
    return(NULL)
  }
  
  # Extract data row
  target_row <- df[target_idx, ]
  sample_names <- colnames(df)[2:ncol(df)]
  values <- as.numeric(target_row[2:ncol(df)])
  
  temp_res <- data.frame(
    Sample = sample_names,
    Value = values,
    Cancer = cancer_name,
    stringsAsFactors = FALSE
  )
  
  # Define Group based on sample suffix
  if (cancer_name %in% cancers_with_normal) {
    # Get the last character of the sample ID
    last_char <- substr(temp_res$Sample, nchar(temp_res$Sample), nchar(temp_res$Sample))
    
    # Logic: 'A' = Normal, 'T' = Tumor (Specific to CPTAC naming convention)
    temp_res$Group <- ifelse(last_char == "A", "Normal", 
                             ifelse(last_char == "T", "Tumor", NA))
    
    # Filter out undefined groups
    temp_res <- temp_res %>% filter(!is.na(Group))
    
  } else {
    # If cancer type doesn't have paired normal in this list, assume all are Tumor
    temp_res$Group <- "Tumor"
  }
  
  return(temp_res)
}

# ==============================================================================
# 3. Process All Cancer Types
# ==============================================================================

for (cancer in cancer_types) {
  res <- extract_rna_data(cancer)
  if (!is.null(res)) {
    plot_data_list[[cancer]] <- res
  }
}

# ==============================================================================
# 4. Data Merging & Validation
# ==============================================================================

if (length(plot_data_list) == 0) stop("No valid data found.")

final_df <- bind_rows(plot_data_list)
final_df <- final_df %>% filter(!is.na(Value))
final_df$Group <- factor(final_df$Group, levels = c("Normal", "Tumor"))

# Update cancer levels based on available data
available_cancers <- unique(final_df$Cancer)
final_levels <- intersect(cancer_types, available_cancers)
final_df$Cancer <- factor(final_df$Cancer, levels = final_levels)

message("Data preparation complete.")

# ==============================================================================
# 4.5 Calculate Signed Meta P-value (Stouffer's Method)
# ==============================================================================
message("Calculating Signed Meta P-value (Stouffer's Method)...")

z_scores <- c()

for (c in final_levels) {
  sub_df <- final_df %>% filter(Cancer == c)
  
  # Calculation requires both Tumor and Normal groups
  if (n_distinct(sub_df$Group) == 2) {
    
    # 1. Wilcoxon test (two-sided)
    res <- wilcox.test(Value ~ Group, data = sub_df)
    p_val <- res$p.value
    
    # Handle extremely small p-values to avoid Inf Z-scores
    if (p_val < 1e-300) p_val <- 1e-300
    
    # 2. Determine direction (Tumor median vs Normal median)
    med_t <- median(sub_df$Value[sub_df$Group == "Tumor"], na.rm=TRUE)
    med_n <- median(sub_df$Value[sub_df$Group == "Normal"], na.rm=TRUE)
    
    # Sign: +1 if Tumor > Normal, -1 if Tumor < Normal
    direction <- sign(med_t - med_n)
    if(direction == 0) direction <- 1 
    
    # 3. Convert to signed Z-score
    z <- qnorm(1 - p_val / 2) * direction
    z_scores <- c(z_scores, z)
  }
}

# 4. Combine Z-scores
if (length(z_scores) > 0) {
  k <- length(z_scores)
  meta_z <- sum(z_scores) / sqrt(k)
  
  # 5. Convert back to P-value
  meta_p_val <- 2 * (1 - pnorm(abs(meta_z)))
  
  # 6. Generate Label
  sign_symbol <- ifelse(meta_z > 0, "(+)", "(-)")
  
  if (meta_p_val < 0.0001) {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " < 0.0001")
  } else {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " = ", sprintf("%.4f", meta_p_val))
  }
} else {
  meta_label <- "Meta P: N/A"
}

message(paste("Calculation Complete:", meta_label))

# ==============================================================================
# 5. Visualization (Style: No black box outlines, White median line)
# ==============================================================================

# Define colors
fill_cols <- c("Tumor" = "#E31A1C", "Normal" = "#377EB8")
line_cols <- c("Tumor" = "#E31A1C", "Normal" = "#377EB8") 

p <- ggplot(final_df, aes(x = Cancer, y = Value, fill = Group, color = Group)) +
  
  # Layer 1: Boxplot Body (fatten = 0 hides the default black median line)
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,
    fatten = 0,             
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.colour = NULL,
    outlier.shape = NA   
  ) +
  
  # Layer 2: Custom White Median Line
  stat_summary(
    fun = median,
    geom = "errorbar",      
    aes(ymax = ..y.., ymin = ..y.., group = Group), 
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,            
    color = "white",        
    size = 0.8              
  ) +
  
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = line_cols) +
  
  theme_classic() +
  
  # Statistical comparison
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", 
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     hide.ns = FALSE,
                     vjust = -0.5,
                     color = "black") + 
  
  # Add Meta P-value annotation
  annotate("text", 
           x = 0.6,                          
           y = max(final_df$Value) * 1.05,  
           label = meta_label, 
           hjust = 0,                        
           size = 5,                         
           fontface = "bold",                
           color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  labs(x = "", y = "TPM (Transcripts Per Million)", 
       title = paste0(target_gene, " mRNA Expression")) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
    panel.border = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save Plot
output_file <- file.path(output_dir, paste0(target_gene, "_Pancancer_mRNA_Boxplot.pdf"))
ggsave(output_file, p, width = 7, height = 6)
message(paste("Plot saved to:", output_file))





















# ==============================================================================
# Script Name: Pancancer_Phosphosite_Analysis.R
# Description: 
#   1. Extracts specific phosphosite abundance from CPTAC Pan-cancer data.
#   2. Matches specific protein/site ID strings (e.g., S156).
#   3. Calculates Signed Meta P-value using Stouffer's Z-score method.
#   4. Generates a boxplot with statistical annotations.
#
# Input Data Requirement:
#   Format: .txt (Tab-separated)
#   Path: ./data/{CancerType}/{CancerType}_phospho_site_abundance_...
#   Structure: First column must contain the specific ID string.
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)       

# ==============================================================================
# 1. Configuration & Preparation
# ==============================================================================

# [Config]: Target Identification
# This string must EXACTLY match the first column in your data files
target_id_string <- "ENSG00000143970.17|ENSP00000391447.3|S156|IPAGKVISPSQKHSK|1"

# [Config]: Labels for Plotting
gene_symbol <- "ASXL2"
site_label <- "S156"

# [Config]: Paths (Relative paths recommended for GitHub)
base_dir <- "./data/CPTACPancancer"
output_dir <- "./results"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# [Config]: Cancer Types
cancer_types <- c("LUAD", "LSCC", "PDAC", "HNSCC", "CCRCC", "UCEC",    
                  "OV", "COAD", "GBM", "BRCA")

# Initialize list to store data
plot_data_list <- list()

message(paste("Target Site ID:", target_id_string))

# ==============================================================================
# 2. Helper Function: Extract Phosphosite Data
# ==============================================================================
extract_phos_data <- function(filepath, cancer_name, group_type) {
  
  # 1. Check if file exists
  if (!file.exists(filepath)) return(NULL)
  
  # 2. Read file using fread (efficient for large phospho files)
  dt <- fread(filepath, header = TRUE, sep = "\t")
  
  # 3. Find target row
  # Assumes column 1 is the Index/ID column
  target_row <- dt[dt[[1]] == target_id_string, ]
  
  # If site not found in this cancer type, return NULL
  if (nrow(target_row) == 0) {
    return(NULL)
  }
  
  # 4. Extract values
  # Assumes data starts from column 2
  sample_names <- colnames(dt)[2:ncol(dt)]
  values <- as.numeric(target_row[1, 2:ncol(dt), with = FALSE])
  
  # 5. Build Data Frame
  df_res <- data.frame(
    Sample = sample_names,
    Value = values,
    Group = group_type,
    Cancer = cancer_name,
    stringsAsFactors = FALSE
  )
  
  # 6. [Critical]: Filter NA values
  # Phosphoproteomics data often contains many NAs; they must be removed.
  df_res <- df_res %>% filter(!is.na(Value))
  
  # If no samples remain after filtering, return NULL
  if (nrow(df_res) == 0) return(NULL)
  
  return(df_res)
}

# ==============================================================================
# 3. Process All Cancer Types
# ==============================================================================

for (cancer in cancer_types) {
  message(paste("Processing:", cancer, "..."))
  dir_path <- file.path(base_dir, cancer)
  
  # Construct file paths
  f_tumor <- file.path(dir_path, paste0(cancer, "_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"))
  f_normal <- file.path(dir_path, paste0(cancer, "_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt"))
  
  dat_t <- extract_phos_data(f_tumor, cancer, "Tumor")
  dat_n <- extract_phos_data(f_normal, cancer, "Normal")
  
  if (!is.null(dat_t)) plot_data_list[[paste0(cancer, "_T")]] <- dat_t
  if (!is.null(dat_n)) plot_data_list[[paste0(cancer, "_N")]] <- dat_n
}

# ==============================================================================
# 4. Data Merging & Validation
# ==============================================================================

if (length(plot_data_list) == 0) {
  stop("No valid data found for this site in any cancer type. Check the ID string or data paths.")
}

final_df <- bind_rows(plot_data_list)

# Double check NA removal
final_df <- final_df %>% filter(!is.na(Value))

# Set Factor Levels
final_df$Group <- factor(final_df$Group, levels = c("Normal", "Tumor"))

# Update cancer levels based on available data
available_cancers <- unique(final_df$Cancer)
final_levels <- intersect(cancer_types, available_cancers)
final_df$Cancer <- factor(final_df$Cancer, levels = final_levels)

message("Data preparation complete. Sample distribution:")
print(table(final_df$Cancer, final_df$Group))

# ==============================================================================
# 4.5 Calculate Signed Meta P-value (Stouffer's Method)
# ==============================================================================
message("Calculating Signed Meta P-value (Stouffer's Method)...")

z_scores <- c()
valid_cancers_for_meta <- 0

for (c in final_levels) {
  sub_df <- final_df %>% filter(Cancer == c)
  
  # Requires both Tumor and Normal samples to calculate statistics
  n_tumor <- sum(sub_df$Group == "Tumor")
  n_normal <- sum(sub_df$Group == "Normal")
  
  if (n_tumor > 0 && n_normal > 0) {
    
    # 1. Wilcoxon Test (two-sided)
    res <- wilcox.test(Value ~ Group, data = sub_df)
    p_val <- res$p.value
    
    if (is.na(p_val)) next 
    if (p_val < 1e-300) p_val <- 1e-300
    
    # 2. Determine Direction
    med_t <- median(sub_df$Value[sub_df$Group == "Tumor"], na.rm=TRUE)
    med_n <- median(sub_df$Value[sub_df$Group == "Normal"], na.rm=TRUE)
    
    direction <- sign(med_t - med_n)
    if(direction == 0) direction <- 1
    
    # 3. Convert to Signed Z-score
    z <- qnorm(1 - p_val / 2) * direction
    z_scores <- c(z_scores, z)
    valid_cancers_for_meta <- valid_cancers_for_meta + 1
  }
}

# 4. Combine Z-scores
if (length(z_scores) > 0) {
  k <- length(z_scores)
  meta_z <- sum(z_scores) / sqrt(k)
  meta_p_val <- 2 * (1 - pnorm(abs(meta_z)))
  
  sign_symbol <- ifelse(meta_z > 0, "(+)", "(-)")
  
  if (meta_p_val < 0.0001) {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " < 0.0001")
  } else {
    meta_label <- paste0("Signed Meta P ", sign_symbol, " = ", sprintf("%.4f", meta_p_val))
  }
} else {
  meta_label <- "Meta P: N/A"
}

message(paste("Calculation Complete:", meta_label))

# ==============================================================================
# 5. Visualization
# ==============================================================================

cols <- c("Tumor" = "#E31A1C", "Normal" = "#377EB8")

p <- ggplot(final_df, aes(x = Cancer, y = Value, fill = Group, color = Group)) +
  
  # Boxplot Body
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,
    fatten = 0,             
    outlier.size = 0.5,
    outlier.alpha = 0.8,
    outlier.colour = NULL,
    outlier.shape = NA   
  ) +
  
  # Custom White Median Line
  stat_summary(
    fun = median,
    geom = "errorbar",      
    aes(ymax = ..y.., ymin = ..y.., group = Group), 
    position = position_dodge2(preserve = "single", padding = 0.2), 
    width = 0.7,            
    color = "white",        
    size = 0.8              
  ) +
  
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  theme_classic() +
  
  # Statistical Comparison
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", 
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     hide.ns = FALSE,
                     vjust = -0.5,
                     color = "black") + 
  
  # Meta P-value Annotation
  annotate("text", 
           x = 0.6,                          
           y = max(final_df$Value) * 1.05,  
           label = meta_label, 
           hjust = 0,                        
           size = 5,                         
           fontface = "bold",                
           color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  labs(x = "", y = "Log2 Intensity", 
       title = paste0(gene_symbol, " (", site_label, ") Phosphorylation")) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
    panel.border = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save Plot
output_filename <- file.path(output_dir, paste0(gene_symbol, "_", site_label, "_Pancancer_Phospho_Boxplot.pdf"))
ggsave(output_filename, p, width = 7, height = 6)

message(paste("Plot saved to:", output_filename))




























# ==============================================================================
# Script Name: Pancancer_Phospho_GSEA.R
# Description: 
#   1. Performs differential analysis (T-test) on Proteomics data based on 
#      Phosphosite expression levels (High vs Low median split).
#   2. Runs GSEA (Gene Set Enrichment Analysis) using Hallmark gene sets.
#   3. Visualizes specific oncogenic and immune pathways in a combined bubble plot.
#
# Input Data Structure:
#   ./data/CPTACPancancer/{CancerType}/{CancerType}_phospho_site_abundance_...
#   ./data/CPTACPancancer/{CancerType}/{CancerType}_proteomics_gene_abundance_...
# ==============================================================================

# ==============================================================================
# 0. Load necessary libraries
# ==============================================================================
library(data.table)
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggplot2)

# ==============================================================================
# 1. Configuration & Preparation
# ==============================================================================

# [Config]: Target Phosphosite Information
# The ID string must match the first column in your data files exactly.
target_phos_id <- "ENSG00000143970.17|ENSP00000391447.3|S156|IPAGKVISPSQKHSK|1"
target_symbol  <- "ASXL2"

# [Config]: Cancer Types to Analyze
cancer_types <- c("LUAD", "LSCC", "PDAC", "HNSCC", "CCRCC", "UCEC", 
                  "OV", "COAD", "GBM", "BRCA")

# [Config]: Paths (Relative paths recommended for GitHub)
base_dir <- "./data/CPTACPancancer"
out_dir  <- "./results/GSEA_Results_Ttest"

# Create output directory if it does not exist
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Prepare Hallmark Gene Sets
message("Loading Hallmark gene sets...")
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_terms <- hallmark_gene_sets[, c("gs_name", "gene_symbol")]
colnames(hallmark_terms) <- c("gs_name", "gene_symbol")

# ==============================================================================
# 2. Core Analysis Function
# ==============================================================================

run_gsea_for_cancer <- function(cancer) {
  message(paste0("\n[", cancer, "] Starting analysis..."))
  dir_path <- file.path(base_dir, cancer)
  
  # --- 2.1 Load Phosphorylation Data & Perform Grouping (Median Split) ---
  f_phos <- file.path(dir_path, paste0(cancer, "_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"))
  if (!file.exists(f_phos)) {
    warning(paste("Phospho file not found for:", cancer))
    return(NULL)
  }
  
  dt_phos <- fread(f_phos, header = TRUE, sep = "\t")
  target_row <- dt_phos[dt_phos[[1]] == target_phos_id, ]
  
  if (nrow(target_row) == 0) {
    message("  - Target phosphosite not found.")
    return(NULL)
  }
  
  phos_cols <- colnames(dt_phos)[2:ncol(dt_phos)]
  phos_vals <- as.numeric(target_row[1, ..phos_cols])
  names(phos_vals) <- phos_cols
  phos_vals <- phos_vals[!is.na(phos_vals)]
  
  if (length(phos_vals) < 10) {
    message("  - Not enough valid samples.")
    return(NULL)
  }
  
  # Grouping: High vs Low based on Median
  median_val <- median(phos_vals)
  high_samples <- names(phos_vals)[phos_vals > median_val]
  low_samples <- names(phos_vals)[phos_vals <= median_val]
  
  message(paste0("  - Grouping: High=", length(high_samples), ", Low=", length(low_samples)))
  
  # --- 2.2 Load Proteomics Data ---
  f_prot <- file.path(dir_path, paste0(cancer, "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))
  if (!file.exists(f_prot)) return(NULL)
  
  dt_prot <- fread(f_prot, header = TRUE, sep = "\t")
  # Clean Ensembl IDs (remove version number)
  dt_prot$ensembl_id <- sub("\\.[0-9]+$", "", dt_prot[[1]])
  
  prot_samples <- colnames(dt_prot)
  valid_high <- intersect(high_samples, prot_samples)
  valid_low <- intersect(low_samples, prot_samples)
  
  if (length(valid_high) < 3 || length(valid_low) < 3) return(NULL)
  
  # --- 2.3 Differential Analysis (T-test) ---
  message("  - Calculating T-statistics for all proteins...")
  
  mat_high <- as.matrix(dt_prot[, ..valid_high])
  mat_low <- as.matrix(dt_prot[, ..valid_low])
  
  # Helper function for robust T-test
  calc_t_stat <- function(x_high, x_low) {
    # Require at least 2 non-NA values per group
    if (sum(!is.na(x_high)) < 2 || sum(!is.na(x_low)) < 2) return(NA)
    tryCatch({
      # Default: Welch t-test
      res <- t.test(x_high, x_low)
      return(res$statistic) # Return t-statistic
    }, error = function(e) return(NA))
  }
  
  # Apply T-test to every row (protein)
  t_stats <- sapply(1:nrow(dt_prot), function(i) {
    calc_t_stat(mat_high[i, ], mat_low[i, ])
  })
  
  res_df <- data.frame(
    ENSEMBL = dt_prot$ensembl_id,
    Score = t_stats, # Using t-statistic for ranking
    stringsAsFactors = FALSE
  )
  res_df <- res_df[!is.na(res_df$Score), ]
  
  # --- 2.4 ID Conversion (Ensembl -> Symbol) ---
  gene_map <- bitr(res_df$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  res_df <- inner_join(res_df, gene_map, by = "ENSEMBL")
  
  # Deduplication: Keep the one with the highest absolute t-statistic
  res_df <- res_df %>% 
    arrange(desc(abs(Score))) %>% 
    distinct(SYMBOL, .keep_all = TRUE)
  
  # --- 2.5 Generate Ranked Gene List ---
  gene_list <- res_df$Score
  names(gene_list) <- res_df$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  message(paste0("  - Genes for GSEA: ", length(gene_list)))
  
  # --- 2.6 Run GSEA ---
  mapped_genes <- intersect(names(gene_list), hallmark_terms$gene_symbol)
  if (length(mapped_genes) < 100) return(NULL)
  
  gsea_res <- tryCatch({
    GSEA(
      geneList = gene_list,
      TERM2GENE = hallmark_terms,
      pvalueCutoff = 1,     # Keep all results for visualization
      pAdjustMethod = "BH",
      verbose = FALSE,
      nPermSimple = 10000,
      seed = 123
    )
  }, error = function(e) return(NULL))
  
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) return(NULL)
  
  # --- 2.7 Save Results ---
  write.csv(gsea_res@result, file.path(out_dir, paste0(cancer, "_GSEA_Ttest_Table.csv")))
  
  # Optional: Save individual dotplot for significant results
  sig_res <- gsea_res@result %>% filter(p.adjust < 0.05)
  if (nrow(sig_res) > 0) {
    p <- dotplot(gsea_res, showCategory = 20, split = ".sign") + facet_grid(.~.sign) +
      ggtitle(paste0(cancer, ": ", target_symbol, " High vs Low (t-stat)")) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(out_dir, paste0(cancer, "_GSEA_Dotplot.pdf")), p, width = 10, height = 8)
  }
  message("  - Completed.")
}

# ==============================================================================
# 3. Execution Loop
# ==============================================================================
for (cancer in cancer_types) {
  run_gsea_for_cancer(cancer)
}

# ==============================================================================
# 4. Visualization: Combined Bubble Plot
#    (Aggregates results from all cancers for specific pathways)
# ==============================================================================

# 4.1 Define Target Pathways and Display Names
pathway_map <- c(
  # --- Pro-tumorigenic / Proliferation ---
  "HALLMARK_E2F_TARGETS"               = "E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"            = "G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1"            = "MYC_V1",
  "HALLMARK_MYC_TARGETS_V2"            = "MYC_V2",
  
  # --- Immune / Inflammation ---
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFN_GAMMA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "IFN_ALPHA",
  "HALLMARK_ALLOGRAFT_REJECTION"       = "ALLOGRAFT",
  "HALLMARK_INFLAMMATORY_RESPONSE"     = "INFLAMMATION"
)

target_pathways <- names(pathway_map)

# 4.2 Read and Merge Data
all_gsea_data <- list()

for (cancer in cancer_types) {
  csv_file <- file.path(out_dir, paste0(cancer, "_GSEA_Ttest_Table.csv"))
  if (file.exists(csv_file)) {
    df <- read.csv(csv_file, row.names = 1, stringsAsFactors = FALSE)
    df_filtered <- df %>% 
      filter(ID %in% target_pathways) %>%
      mutate(Cancer = cancer) 
    if (nrow(df_filtered) > 0) all_gsea_data[[cancer]] <- df_filtered
  }
}

if (length(all_gsea_data) == 0) stop("No data found for target pathways.")
plot_df <- bind_rows(all_gsea_data)

# 4.3 Data Preprocessing for Plotting
plot_df <- plot_df %>%
  mutate(
    CoreCount = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = CoreCount / setSize,
    logP = -log10(ifelse(p.adjust < 1e-10, 1e-10, p.adjust)),
    ShortName = pathway_map[ID],
    ShortName = factor(ShortName, levels = c(
      "E2F_TARGETS", "G2M_CHECKPOINT", "MYC_V1","MYC_V2",
      "IFN_GAMMA", "IFN_ALPHA", "ALLOGRAFT","INFLAMMATION"
    ))
  )

# Sorting Cancers based on E2F_TARGETS GeneRatio
e2f_data <- plot_df %>%
  filter(ShortName == "E2F_TARGETS") %>%
  select(Cancer, GeneRatio) %>%
  arrange(GeneRatio) 

ordered_cancers <- e2f_data$Cancer
all_cancers <- unique(plot_df$Cancer)
missing_cancers <- setdiff(all_cancers, ordered_cancers)
final_levels <- c(missing_cancers, ordered_cancers)

plot_df$Cancer <- factor(plot_df$Cancer, levels = final_levels)

# ==============================================================================
# 5. Generate Bubble Plot (ggplot2)
# ==============================================================================

max_nes <- max(abs(plot_df$NES), na.rm = TRUE)

p_bubble <- ggplot(plot_df, aes(x = GeneRatio, y = Cancer)) +
  
  geom_point(aes(size = logP, fill = NES), shape = 21, color = "black", stroke = 0.6) +
  
  facet_grid(. ~ ShortName, scales = "free_x") + 
  
  scale_fill_gradient2(
    low = "#377EB8", mid = "white", high = "#E31A1C", 
    midpoint = 0, limits = c(-max_nes, max_nes), name = "NES"
  ) +
  
  scale_size_continuous(range = c(3, 10), name = "-log10(adj.P)") +
  
  # [Customization]: Force x-axis ticks to be 0.1, 0.2, ...
  scale_x_continuous(
    breaks = seq(0, 1, 0.1),                
    expand = expansion(mult = c(0.1, 0.15)) 
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"), 
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    
    panel.spacing.x = unit(0.3, "cm"), 
    
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "GeneRatio", 
    y = "",
    title = paste0(target_symbol, " Oncogenic & Immune Pathways")
  )

# ==============================================================================
# 6. Save Plot
# ==============================================================================
output_file <- file.path(out_dir, "ASXL2_S156_Combined_Bubble.pdf")
ggsave(output_file, p_bubble, width = 18, height = 6.5)

message(paste("Combined plot saved to:", output_file))
