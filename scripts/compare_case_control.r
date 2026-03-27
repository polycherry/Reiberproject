# Case vs Control Comparison Analysis
# Compare baseline measurements between Case and Control groups

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

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
# IgG Index = (CSF IgG / Serum IgG) / (CSF Albumin / Serum Albumin)
# Also known as IgG Quotient Ratio or QIgG/QAlb
df$IgG_Index <- (df$CSF_IgG / df$Serum_IgG) / (df$CSF_Albumin / df$Serum_Albumin)

# Albumin Quotient (QAlb)
df$Albumin_Index <- df$CSF_Albumin / df$Serum_Albumin

# IgG Quotient (QIgG)
df$IgG_Quotient <- df$CSF_IgG / df$Serum_IgG

# Filter to baseline only
df_baseline <- df[df$Type == "Baseline" & !is.na(df$Group), ]

# Remove rows with missing values for comparison
df_baseline <- df_baseline[complete.cases(df_baseline[, c("CSF_IgG", "Serum_IgG", "CSF_Albumin", "Serum_Albumin")]), ]

# Print summary statistics
cat("=== BASELINE DATA SUMMARY ===\n")
cat(sprintf("Cases: %d, Controls: %d\n\n", 
            sum(df_baseline$Group == "Case"), 
            sum(df_baseline$Group == "Control")))

# Function to compare groups and create plots
compare_groups <- function(data, variable_name, variable_label, log_scale = FALSE) {
  
  # Extract data for cases and controls
  case_data <- data[data$Group == "Case", variable_name]
  control_data <- data[data$Group == "Control", variable_name]
  
  # Remove NA values
  case_data <- case_data[!is.na(case_data)]
  control_data <- control_data[!is.na(control_data)]
  
  # Statistical test (Wilcoxon rank-sum test / Mann-Whitney U test)
  if (length(case_data) > 0 && length(control_data) > 0) {
    test_result <- wilcox.test(case_data, control_data)
    
    # Print statistics
    cat(sprintf("\n--- %s ---\n", variable_label))
    cat(sprintf("Case: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n", 
                length(case_data), median(case_data), 
                quantile(case_data, 0.25), quantile(case_data, 0.75)))
    cat(sprintf("Control: n=%d, median=%.3f, IQR=[%.3f, %.3f]\n", 
                length(control_data), median(control_data),
                quantile(control_data, 0.25), quantile(control_data, 0.75)))
    cat(sprintf("Wilcoxon test p-value: %.4f %s\n", 
                test_result$p.value,
                ifelse(test_result$p.value < 0.05, "***", "")))
    
    # Create plot data
    plot_data <- data.frame(
      Group = data$Group,
      Value = data[[variable_name]]
    )
    plot_data <- plot_data[!is.na(plot_data$Value), ]
    
    # Create boxplot with jittered points
    p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 2.5, alpha = 0.7) +
      scale_fill_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      labs(
        title = variable_label,
        subtitle = sprintf("p = %.4f (Wilcoxon test)", test_result$p.value),
        x = "Group",
        y = variable_label
      )
    
    # Add log scale if requested
    if (log_scale) {
      p <- p + scale_y_log10()
    }
    
    # Save plot
    filename <- file.path("results", paste0("plot_", gsub(" ", "_", variable_name), ".png"))
    ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
    cat(sprintf("Plot saved: %s\n", filename))
    
    # Return results for composite plot
    return(list(
      variable = variable_label,
      p_value = test_result$p.value,
      plot = p
    ))
  }
  return(NULL)
}

# Compare all variables
results <- list()
results[[1]] <- compare_groups(df_baseline, "CSF_Protein", "CSF Protein (mg/L)")
results[[2]] <- compare_groups(df_baseline, "CSF_IgG", "CSF IgG (mg/L)")
results[[3]] <- compare_groups(df_baseline, "Serum_IgG", "Serum IgG (g/L)")
results[[4]] <- compare_groups(df_baseline, "CSF_Albumin", "CSF Albumin (mg/L)")
results[[5]] <- compare_groups(df_baseline, "Serum_Albumin", "Serum Albumin (g/L)")
results[[6]] <- compare_groups(df_baseline, "Age", "Age (years)")

# Compare calculated indices
results[[7]] <- compare_groups(df_baseline, "IgG_Index", "IgG Index (QIgG/QAlb)", log_scale = TRUE)
results[[8]] <- compare_groups(df_baseline, "Albumin_Index", "Albumin Quotient (QAlb)", log_scale = TRUE)
results[[9]] <- compare_groups(df_baseline, "IgG_Quotient", "IgG Quotient (QIgG)", log_scale = TRUE)

# Remove NULL results
results <- results[!sapply(results, is.null)]

# Create composite plot showing all comparisons
cat("\n=== Creating composite plot ===\n")

# Prepare data for composite plot
composite_data <- data.frame(
  Variable = sapply(results, function(x) x$variable),
  P_Value = sapply(results, function(x) x$p_value),
  Significant = sapply(results, function(x) x$p_value < 0.05)
)

# Order by p-value
composite_data <- composite_data[order(composite_data$P_Value), ]
composite_data$Variable <- factor(composite_data$Variable, levels = composite_data$Variable)

# Create composite plot
p_composite <- ggplot(composite_data, aes(x = Variable, y = -log10(P_Value), fill = Significant)) +
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
    title = "Case vs Control Comparisons - Baseline",
    subtitle = "Dashed line indicates p = 0.05 significance threshold",
    x = "Variable",
    y = "-log10(p-value)",
    fill = ""
  )

# Save composite plot
ggsave("results/plot_composite_all_comparisons.png", plot = p_composite, 
       width = 10, height = 6, dpi = 300)
cat("Composite plot saved: results/plot_composite_all_comparisons.png\n")

# --- Create Summary Table ---
cat("\n=== Creating Summary Table ===\n")

# Variables to summarize
variables <- c("CSF_Protein", "CSF_IgG", "Serum_IgG", "CSF_Albumin", 
               "Serum_Albumin", "Age", "IgG_Index", "Albumin_Index", "IgG_Quotient")
variable_labels <- c("CSF Protein (mg/L)", "CSF IgG (mg/L)", "Serum IgG (g/L)", 
                     "CSF Albumin (mg/L)", "Serum Albumin (g/L)", "Age (years)",
                     "IgG Index", "Albumin Quotient (QAlb)", "IgG Quotient (QIgG)")

# Create summary table
summary_table <- data.frame(
  Parameter = character(),
  Case_Mean_SD = character(),
  Control_Mean_SD = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(variables)) {
  var <- variables[i]
  
  # Extract data
  case_data <- df_baseline[df_baseline$Group == "Case", var]
  control_data <- df_baseline[df_baseline$Group == "Control", var]
  
  # Remove NA
  case_data <- case_data[!is.na(case_data)]
  control_data <- control_data[!is.na(control_data)]
  
  if (length(case_data) > 0 && length(control_data) > 0) {
    # Calculate mean and SD
    case_mean <- mean(case_data)
    case_sd <- sd(case_data)
    control_mean <- mean(control_data)
    control_sd <- sd(control_data)
    
    # Wilcoxon test
    test_result <- wilcox.test(case_data, control_data)
    
    # Add to table
    summary_table <- rbind(summary_table, data.frame(
      Parameter = variable_labels[i],
      Case_Mean_SD = sprintf("%.3f (%.3f)", case_mean, case_sd),
      Control_Mean_SD = sprintf("%.3f (%.3f)", control_mean, control_sd),
      P_Value = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# Add significance stars
summary_table$Significance <- ifelse(summary_table$P_Value < 0.001, "***",
                                    ifelse(summary_table$P_Value < 0.01, "**",
                                          ifelse(summary_table$P_Value < 0.05, "*", "")))

# Format p-value column
summary_table$P_Value_Formatted <- sprintf("%.4f", summary_table$P_Value)

# Final table for display
final_table <- summary_table[, c("Parameter", "Case_Mean_SD", "Control_Mean_SD", "P_Value_Formatted", "Significance")]
colnames(final_table) <- c("Parameter", "Case Mean (SD)", "Control Mean (SD)", "P-Value", "Sig.")

# Print table
cat("\n=== CASE vs CONTROL COMPARISON (BASELINE ONLY) ===\n")
cat(sprintf("Case n = %d, Control n = %d\n\n", 
            sum(df_baseline$Group == "Case"), 
            sum(df_baseline$Group == "Control")))
print(final_table, row.names = FALSE)
cat("\nStatistical test: Wilcoxon rank-sum test (Mann-Whitney U)\n")
cat("* p < 0.05, ** p < 0.01, *** p < 0.001\n")

# Save to CSV
write.csv(final_table, "results/case_control_comparison_summary.csv", row.names = FALSE)
cat("\nSummary table saved to: results/case_control_comparison_summary.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
