# ---------------------------
# Script: Correlate CSF Parameters with Hemoglobin OD
# ---------------------------
# Purpose: Analyze correlations between CSF parameters (protein, IgG, albumin, etc.)
#          and hemoglobin optical density (OD) measurements
# ---------------------------

# Load required libraries
library(readxl)
library(ggplot2)
library(corrplot)
library(Hmisc)  # For correlation with p-values
library(tidyr)
library(dplyr)

# --- Load Hemoglobin OD Data ---
hb_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_07_CSF_Hb_ELISA/IRIS samples Hb OD.xlsx"
hb_data <- read_excel(hb_path)

cat("Hemoglobin OD data loaded:\n")
cat(sprintf("  Dimensions: %d rows, %d columns\n", nrow(hb_data), ncol(hb_data)))
cat(sprintf("  Columns: %s\n", paste(colnames(hb_data), collapse = ", ")))

# --- Load CSF Data (using the same source as run_reiberplot_ondata.r) ---
data_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/Merged_CSF_Serum_Table_CORRECTED.xlsx"
raw_data <- read_excel(data_path)

# Load age metadata
age_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/metadata_subjectID_age.csv"
age_data <- read.csv(age_path)

# Prepare CSF data frame
csf_data <- data.frame(
  patient_id = raw_data$Sample,
  group = raw_data$Group,
  type = raw_data$Type,
  csf_igg = as.numeric(raw_data$CSF_IgG) / 1000,  # Convert to g/L
  s_igg = as.numeric(raw_data$Serum_IgG),
  csf_alb = as.numeric(raw_data$CSF_Albumin) / 1000,  # Convert to g/L
  s_alb = as.numeric(raw_data$Serum_Albumin)
)

# Merge age data
colnames(age_data) <- c("patient_id", "age")
csf_data <- merge(csf_data, age_data, by = "patient_id", all.x = TRUE)

# Calculate quotients
csf_data$QAlb <- csf_data$csf_alb / csf_data$s_alb
csf_data$QIgG <- csf_data$csf_igg / csf_data$s_igg

cat("\nCSF data loaded:\n")
cat(sprintf("  Dimensions: %d rows, %d columns\n", nrow(csf_data), ncol(csf_data)))

# --- Prepare and merge data ---
# Rename Sample...4 column in hb_data to patient_id for merging
hb_data_clean <- hb_data
colnames(hb_data_clean)[colnames(hb_data_clean) == "Sample...4"] <- "patient_id"

# Keep only relevant columns from Hb data
hb_data_clean <- hb_data_clean[, c("patient_id", "OD")]

# Merge CSF and Hb data by patient_id
merged_data <- merge(csf_data, hb_data_clean, by = "patient_id", all = FALSE)

cat(sprintf("\nMerged data: %d samples with both CSF and Hb measurements\n", nrow(merged_data)))

# Check if merge was successful
if (nrow(merged_data) == 0) {
  stop("No matching samples found between CSF and Hb data. Check sample ID formats.")
}

# --- Select parameters for correlation analysis ---
# CSF parameters of interest
csf_params <- c("csf_igg", "s_igg", "csf_alb", "s_alb", "QAlb", "QIgG", "age")

# Hemoglobin OD parameter
hb_params <- c("OD")

cat(sprintf("\nCSF parameters: %s\n", paste(csf_params, collapse = ", ")))
cat(sprintf("Hemoglobin parameter: %s\n", paste(hb_params, collapse = ", ")))

# Create subset with complete cases
analysis_params <- c(csf_params, hb_params)
analysis_data <- merged_data[, c("patient_id", "group", "type", analysis_params)]

# Remove rows with any missing values in the parameters
analysis_data_complete <- analysis_data[complete.cases(analysis_data[, analysis_params]), ]

cat(sprintf("\nComplete cases for analysis: %d samples\n", nrow(analysis_data_complete)))

# --- Calculate Correlations ---patient_id", "group", "t
cor_data <- analysis_data_complete[, analysis_params]

# Convert all columns to numeric
cor_data <- as.data.frame(lapply(cor_data, as.numeric))

# Calculate correlation matrix with p-values
cor_result <- rcorr(as.matrix(cor_data), type = "pearson")
cor_matrix <- cor_result$r
p_matrix <- cor_result$P

# Extract correlations between CSF parameters and Hb parameters
csf_hb_correlations <- cor_matrix[csf_params, hb_params, drop = FALSE]
csf_hb_pvalues <- p_matrix[csf_params, hb_params, drop = FALSE]

# --- Create correlation summary table ---
cor_summary <- data.frame()
for (csf in csf_params) {
  for (hb in hb_params) {
    if (!is.na(csf_hb_correlations[csf, hb])) {
      cor_summary <- rbind(cor_summary, data.frame(
        CSF_Parameter = csf,
        Hb_Parameter = hb,
        Correlation = csf_hb_correlations[csf, hb],
        P_value = csf_hb_pvalues[csf, hb],
        Significant = ifelse(csf_hb_pvalues[csf, hb] < 0.05, "Yes", "No")
      ))
    }
  }
}

# Sort by absolute correlation strength
cor_summary <- cor_summary[order(abs(cor_summary$Correlation), decreasing = TRUE), ]

# Print correlation summary
cat("\n=== Correlation Summary (CSF Parameters vs Hemoglobin OD) ===\n")
print(cor_summary, row.names = FALSE)

# Save correlation summary
write.csv(cor_summary, "results/csf_hb_correlation_summary.csv", row.names = FALSE)
cat("\nCorrelation summary saved to results/csf_hb_correlation_summary.csv\n")

# --- Visualizations ---

# Correlation heatmap (full correlation matrix)
pdf("results/csf_hb_correlation_heatmap.pdf", width = 10, height = 8)
corrplot(cor_matrix, 
         method = "color", 
         type = "upper",
         tl.col = "black", 
         tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         addCoef.col = "black",
         number.cex = 0.6,
         cl.cex = 0.8,
         tl.cex = 0.8,
         title = "Correlation Matrix: CSF Parameters and Hemoglobin OD",
         mar = c(0,0,2,0))
dev.off()
cat("Correlation heatmap saved to results/csf_hb_correlation_heatmap.pdf\n")

# --- Function to plot OD vs any CSF parameter ---
plot_od_vs_parameter <- function(parameter, data = analysis_data_complete, 
                                 save_plot = TRUE, filename = NULL) {
  # Check if parameter exists
  if (!parameter %in% colnames(data)) {
    stop(sprintf("Parameter '%s' not found in data. Available parameters: %s", 
                 parameter, paste(colnames(data), collapse = ", ")))
  }
  
  # Get data
  x <- as.numeric(data$OD)
  y <- as.numeric(data[[parameter]])
  
  # Remove NA values
  valid_idx <- !is.na(x) & !is.na(y)
  x_clean <- x[valid_idx]
  y_clean <- y[valid_idx]
  group_clean <- data$group[valid_idx]
  
  # Calculate correlation
  cor_test <- cor.test(x_clean, y_clean, method = "pearson")
  
  # Create ggplot
  p <- ggplot(data[valid_idx, ], aes(x = OD, y = .data[[parameter]], color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
    labs(title = sprintf("Hemoglobin OD vs %s", parameter),
         subtitle = sprintf("Pearson r = %.3f, p = %.4f, n = %d", 
                          cor_test$estimate, cor_test$p.value, length(x_clean)),
         x = "Hemoglobin OD",
         y = parameter,
         color = "Group") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12))
  
  # Save plot if requested
  if (save_plot) {
    if (is.null(filename)) {
      filename <- sprintf("results/od_vs_%s.pdf", gsub("[^A-Za-z0-9]", "_", parameter))
    }
    ggsave(filename, plot = p, width = 8, height = 6)
    cat(sprintf("Plot saved to %s\n", filename))
  }
  
  return(p)
}

cat("\n=== Analysis Complete ===\n")
cat("All results saved to results/ directory\n")

plot_od_vs_parameter("QAlb")
plot_od_vs_parameter("csf_igg")
