# Paired Baseline vs Followup Comparison
# Compare measurements between baseline and followup for patients with both timepoints

library(readxl)
library(ggplot2)
library(dplyr)

# Load actual data from Excel file
data_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/Merged_CSF_Serum_Table_CORRECTED.xlsx"
raw_data <- read_excel(data_path)

# Load age metadata
age_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/metadata_subjectID_age.csv"
age_data <- read.csv(age_path)

# Prepare data frame
df <- data.frame(
  Sample = raw_data$Sample,
  Group = raw_data$Group,
  Type = raw_data$Type,
  CSF_Protein = as.numeric(raw_data$CSF_Protein),
  CSF_IgG = as.numeric(raw_data$CSF_IgG),
  Serum_IgG = as.numeric(raw_data$Serum_IgG),
  CSF_Albumin = as.numeric(raw_data$CSF_Albumin),
  Serum_Albumin = as.numeric(raw_data$Serum_Albumin)
)

# Merge age data
colnames(age_data) <- c("Sample", "Age")
df <- merge(df, age_data, by = "Sample", all.x = TRUE)

# Calculate indices
df$IgG_Index <- (df$CSF_IgG / df$Serum_IgG) / (df$CSF_Albumin / df$Serum_Albumin)
df$Albumin_Index <- df$CSF_Albumin / df$Serum_Albumin
df$IgG_Quotient <- df$CSF_IgG / df$Serum_IgG

# Extract patient ID (remove letters like 'a' or 'b' suffix for pairing)
df$PatientID <- gsub("[a-z]$", "", df$Sample)

# Find patients who have both baseline and followup
baseline_patients <- unique(df$PatientID[df$Type == "Baseline"])
followup_patients <- unique(df$PatientID[df$Type == "Followup"])
paired_patients <- intersect(baseline_patients, followup_patients)

cat(sprintf("Found %d patients with both baseline and followup measurements\n", length(paired_patients)))

# Create output directory for paired comparisons
output_dir <- "results/baseline_vs_followup"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Filter to only paired patients
df_paired <- df[df$PatientID %in% paired_patients, ]

# Separate baseline and followup
df_baseline <- df_paired[df_paired$Type == "Baseline", ]
df_followup <- df_paired[df_paired$Type == "Followup", ]

# Merge baseline and followup data
df_merged <- merge(df_baseline, df_followup, 
                   by = "PatientID", 
                   suffixes = c("_Baseline", "_Followup"))

# Remove rows with missing values
df_merged <- df_merged[complete.cases(df_merged[, c("CSF_IgG_Baseline", "CSF_IgG_Followup")]), ]

cat(sprintf("Total paired samples for analysis: %d\n\n", nrow(df_merged)))

# Function to create paired before-after plot with connecting lines
create_paired_plot <- function(data, var_baseline, var_followup, variable_label, log_scale = FALSE, output_dir = "results") {
  
  # Prepare data for plotting
  plot_data <- data.frame(
    PatientID = rep(data$PatientID, 2),
    Group = rep(data$Group_Baseline, 2),
    Timepoint = rep(c("Baseline", "Followup"), each = nrow(data)),
    Value = c(data[[var_baseline]], data[[var_followup]])
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$Value), ]
  
  # Make sure timepoint is a factor with correct order
  plot_data$Timepoint <- factor(plot_data$Timepoint, levels = c("Baseline", "Followup"))
  
  # Statistical test (paired Wilcoxon signed-rank test)
  baseline_vals <- data[[var_baseline]]
  followup_vals <- data[[var_followup]]
  
  # Remove pairs with NA
  valid_pairs <- !is.na(baseline_vals) & !is.na(followup_vals)
  baseline_vals <- baseline_vals[valid_pairs]
  followup_vals <- followup_vals[valid_pairs]
  
  if (length(baseline_vals) > 0) {
    test_result <- wilcox.test(baseline_vals, followup_vals, paired = TRUE)
    
    # Print statistics
    cat(sprintf("\n--- %s ---\n", variable_label))
    cat(sprintf("Baseline: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n", 
                length(baseline_vals), median(baseline_vals),
                quantile(baseline_vals, 0.25), quantile(baseline_vals, 0.75)))
    cat(sprintf("Followup: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n", 
                length(followup_vals), median(followup_vals),
                quantile(followup_vals, 0.25), quantile(followup_vals, 0.75)))
    cat(sprintf("Paired Wilcoxon test p-value: %.4f %s\n", 
                test_result$p.value,
                ifelse(test_result$p.value < 0.05, "***", "")))
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Timepoint, y = Value)) +
      # Add connecting lines for each patient
      geom_line(aes(group = PatientID, color = Group), alpha = 0.5, linewidth = 0.8) +
      # Add points
      geom_point(aes(color = Group), size = 3, alpha = 0.8) +
      # Add boxplot overlay
      geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.3, fill = NA, color = "black") +
      scale_color_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "top"
      ) +
      labs(
        title = variable_label,
        subtitle = sprintf("Paired comparison (n=%d), p = %.4f", 
                          length(baseline_vals), test_result$p.value),
        x = "Timepoint",
        y = variable_label,
        color = "Group"
      )
    
    # Add log scale if requested
    if (log_scale) {
      p <- p + scale_y_log10()
    }
    
    # Save plot
    var_name <- gsub("/", "", variable_label)  # Remove slashes first
    var_name <- gsub("[()]", "", var_name)     # Remove parentheses
    var_name <- gsub(" ", "_", var_name)       # Replace spaces with underscores
    var_name <- gsub("__+", "_", var_name)     # Remove multiple underscores
    filename <- file.path(output_dir, paste0("paired_", var_name, ".png"))
    
    ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
    cat(sprintf("Plot saved: %s\n", filename))
    
    return(list(
      variable = variable_label,
      p_value = test_result$p.value,
      n_pairs = length(baseline_vals)
    ))
  }
  return(NULL)
}

# Create paired plots for all variables
cat("=== PAIRED BASELINE vs FOLLOWUP ANALYSIS ===\n")

results <- list()
results[[1]] <- create_paired_plot(df_merged, "CSF_Protein_Baseline", "CSF_Protein_Followup", 
                                    "CSF Protein (mg/L)", output_dir = output_dir)
results[[2]] <- create_paired_plot(df_merged, "CSF_IgG_Baseline", "CSF_IgG_Followup", 
                                    "CSF IgG (mg/L)", output_dir = output_dir)
results[[3]] <- create_paired_plot(df_merged, "Serum_IgG_Baseline", "Serum_IgG_Followup", 
                                    "Serum IgG (g/L)", output_dir = output_dir)
results[[4]] <- create_paired_plot(df_merged, "CSF_Albumin_Baseline", "CSF_Albumin_Followup", 
                                    "CSF Albumin (mg/L)", output_dir = output_dir)
results[[5]] <- create_paired_plot(df_merged, "Serum_Albumin_Baseline", "Serum_Albumin_Followup", 
                                    "Serum Albumin (g/L)", output_dir = output_dir)
results[[6]] <- create_paired_plot(df_merged, "IgG_Index_Baseline", "IgG_Index_Followup", 
                                    "IgG Index", log_scale = TRUE, output_dir = output_dir)
results[[7]] <- create_paired_plot(df_merged, "Albumin_Index_Baseline", "Albumin_Index_Followup", 
                                    "Albumin Quotient (QAlb)", log_scale = TRUE, output_dir = output_dir)
results[[8]] <- create_paired_plot(df_merged, "IgG_Quotient_Baseline", "IgG_Quotient_Followup", 
                                    "IgG Quotient (QIgG)", log_scale = TRUE, output_dir = output_dir)

# Remove NULL results
results <- results[!sapply(results, is.null)]

# Create summary plot
cat("\n=== Creating summary plot ===\n")

summary_data <- data.frame(
  Variable = sapply(results, function(x) x$variable),
  P_Value = sapply(results, function(x) x$p_value),
  N_Pairs = sapply(results, function(x) x$n_pairs),
  Significant = sapply(results, function(x) x$p_value < 0.05)
)

# Order by p-value
summary_data <- summary_data[order(summary_data$P_Value), ]
summary_data$Variable <- factor(summary_data$Variable, levels = summary_data$Variable)

# Create summary plot
p_summary <- ggplot(summary_data, aes(x = Variable, y = -log10(P_Value), fill = Significant)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "#D55E00"),
                    labels = c("Not Significant", "Significant (p < 0.05)")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  ) +
  labs(
    title = "Baseline vs Followup Paired Comparisons",
    subtitle = "Dashed line indicates p = 0.05 significance threshold",
    x = "Variable",
    y = "-log10(p-value)",
    fill = ""
  )

ggsave(file.path(output_dir, "paired_summary_all_comparisons.png"), plot = p_summary, 
       width = 10, height = 6, dpi = 300)
cat(sprintf("Summary plot saved: %s/paired_summary_all_comparisons.png\n", output_dir))

cat("\n=== ANALYSIS COMPLETE ===\n")
