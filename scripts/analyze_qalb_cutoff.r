# Individual QAlb Analysis - Compare to Age-Adjusted Cutoffs
# Assess blood-CSF barrier dysfunction relative to individual patient age

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
  CSF_Albumin = as.numeric(raw_data$CSF_Albumin),
  Serum_Albumin = as.numeric(raw_data$Serum_Albumin)
)

# Merge age data
colnames(age_data) <- c("Sample", "Age")
df <- merge(df, age_data, by = "Sample", all.x = TRUE)

# Calculate QAlb for each patient
df$QAlb <- (df$CSF_Albumin / 1000) / df$Serum_Albumin  # Convert CSF_Albumin to g/L

# Calculate individualized QAlb cutoff based on age
# Formula from Reibergram: qalb_cutoff = (4 + age/15) * 0.001
df$QAlb_Cutoff <- (4 + df$Age / 15) * 0.001

# Calculate ratio: QAlb / Cutoff
# Ratio > 1 indicates blood-CSF barrier dysfunction
df$QAlb_Ratio <- df$QAlb / df$QAlb_Cutoff

# Calculate absolute distance from cutoff
df$QAlb_Distance <- df$QAlb - df$QAlb_Cutoff

# Remove rows with missing values
df_complete <- df[complete.cases(df[, c("QAlb", "QAlb_Cutoff", "Group", "Type")]), ]

# Filter to baseline only
df_baseline <- df_complete[df_complete$Type == "Baseline", ]

cat("=== INDIVIDUALIZED QAlb ANALYSIS ===\n")
cat(sprintf("Total baseline samples: %d\n", nrow(df_baseline)))
cat(sprintf("Cases: %d, Controls: %d\n\n", 
            sum(df_baseline$Group == "Case"), 
            sum(df_baseline$Group == "Control")))

# Statistical comparison
case_ratio <- df_baseline$QAlb_Ratio[df_baseline$Group == "Case"]
control_ratio <- df_baseline$QAlb_Ratio[df_baseline$Group == "Control"]

test_result <- wilcox.test(case_ratio, control_ratio)

cat("--- QAlb/Cutoff Ratio ---\n")
cat(sprintf("Cases: median=%.3f, IQR=[%.3f, %.3f]\n", 
            median(case_ratio, na.rm = TRUE),
            quantile(case_ratio, 0.25, na.rm = TRUE),
            quantile(case_ratio, 0.75, na.rm = TRUE)))
cat(sprintf("Controls: median=%.3f, IQR=[%.3f, %.3f]\n", 
            median(control_ratio, na.rm = TRUE),
            quantile(control_ratio, 0.25, na.rm = TRUE),
            quantile(control_ratio, 0.75, na.rm = TRUE)))
cat(sprintf("Wilcoxon test p-value: %.4f %s\n\n", 
            test_result$p.value,
            ifelse(test_result$p.value < 0.05, "***", "")))

# Count patients above cutoff
cat("--- Barrier Dysfunction (QAlb > Cutoff) ---\n")
cat(sprintf("Cases: %d/%d (%.1f%%)\n", 
            sum(df_baseline$QAlb_Ratio[df_baseline$Group == "Case"] > 1),
            sum(df_baseline$Group == "Case"),
            100 * sum(df_baseline$QAlb_Ratio[df_baseline$Group == "Case"] > 1) / sum(df_baseline$Group == "Case")))
cat(sprintf("Controls: %d/%d (%.1f%%)\n\n", 
            sum(df_baseline$QAlb_Ratio[df_baseline$Group == "Control"] > 1),
            sum(df_baseline$Group == "Control"),
            100 * sum(df_baseline$QAlb_Ratio[df_baseline$Group == "Control"] > 1) / sum(df_baseline$Group == "Control")))

# Create output directory
output_dir <- "results/qalb_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plot 1: QAlb Ratio boxplot with jittered points
p1 <- ggplot(df_baseline, aes(x = Group, y = QAlb_Ratio, fill = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "QAlb Relative to Age-Adjusted Cutoff",
    subtitle = sprintf("p = %.4f (Wilcoxon test), Red line = cutoff", test_result$p.value),
    x = "Group",
    y = "QAlb / Age-Adjusted Cutoff"
  ) +
  annotate("text", x = 1.5, y = max(df_baseline$QAlb_Ratio, na.rm = TRUE) * 0.95,
           label = "Ratio > 1 = Barrier dysfunction", color = "red", fontface = "bold")

ggsave(file.path(output_dir, "qalb_ratio_boxplot.png"), plot = p1, 
       width = 8, height = 6, dpi = 300)
cat("Plot saved: qalb_ratio_boxplot.png\n")

# Plot 2: Individual patient plot showing QAlb vs Age-adjusted cutoff
p2 <- ggplot(df_baseline, aes(x = Age, y = QAlb, color = Group)) +
  # Add cutoff line
  geom_line(aes(y = QAlb_Cutoff), color = "red", linewidth = 1.5, linetype = "dashed") +
  # Add points
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  ) +
  labs(
    title = "QAlb vs Age-Adjusted Cutoff",
    subtitle = "Red dashed line shows age-dependent upper limit of normal",
    x = "Age (years)",
    y = "QAlb (CSF Albumin / Serum Albumin)",
    color = "Group"
  )

ggsave(file.path(output_dir, "qalb_vs_age_cutoff.png"), plot = p2, 
       width = 10, height = 6, dpi = 300)
cat("Plot saved: qalb_vs_age_cutoff.png\n")

# Plot 3: Distance from cutoff (absolute difference)
p3 <- ggplot(df_baseline, aes(x = Group, y = QAlb_Distance * 1000, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "Distance from Age-Adjusted QAlb Cutoff",
    subtitle = "Positive values indicate barrier dysfunction",
    x = "Group",
    y = "QAlb - Cutoff (×10⁻³)"
  )

ggsave(file.path(output_dir, "qalb_distance_from_cutoff.png"), plot = p3, 
       width = 8, height = 6, dpi = 300)
cat("Plot saved: qalb_distance_from_cutoff.png\n")

# Plot 4: Scatter plot showing relationship between QAlb ratio and age
p4 <- ggplot(df_baseline, aes(x = Age, y = QAlb_Ratio, color = Group, shape = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
  scale_shape_manual(values = c("Case" = 16, "Control" = 17)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  ) +
  labs(
    title = "QAlb Ratio vs Age",
    subtitle = "Normalized to individual age-adjusted cutoffs",
    x = "Age (years)",
    y = "QAlb / Age-Adjusted Cutoff",
    color = "Group",
    shape = "Group"
  ) +
  annotate("text", x = mean(df_baseline$Age, na.rm = TRUE), 
           y = max(df_baseline$QAlb_Ratio, na.rm = TRUE) * 0.95,
           label = "Above red line = Barrier dysfunction", 
           color = "red", fontface = "bold")

ggsave(file.path(output_dir, "qalb_ratio_vs_age.png"), plot = p4, 
       width = 10, height = 6, dpi = 300)
cat("Plot saved: qalb_ratio_vs_age.png\n")

# Plot 5: Violin plot showing distribution
p5 <- ggplot(df_baseline, aes(x = Group, y = QAlb_Ratio, fill = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.5) +
  scale_fill_manual(values = c("Case" = "#D55E00", "Control" = "#56B4E9")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "Distribution of QAlb Ratios",
    subtitle = sprintf("Cases vs Controls (p = %.4f)", test_result$p.value),
    x = "Group",
    y = "QAlb / Age-Adjusted Cutoff"
  )

ggsave(file.path(output_dir, "qalb_ratio_violin.png"), plot = p5, 
       width = 8, height = 6, dpi = 300)
cat("Plot saved: qalb_ratio_violin.png\n")

# Create a summary table
summary_table <- df_baseline %>%
  group_by(Group) %>%
  summarise(
    N = n(),
    Mean_Age = mean(Age, na.rm = TRUE),
    Median_QAlb = median(QAlb, na.rm = TRUE),
    Median_Cutoff = median(QAlb_Cutoff, na.rm = TRUE),
    Median_Ratio = median(QAlb_Ratio, na.rm = TRUE),
    N_Above_Cutoff = sum(QAlb_Ratio > 1),
    Percent_Above_Cutoff = 100 * sum(QAlb_Ratio > 1) / n()
  )

cat("\n=== SUMMARY TABLE ===\n")
print(summary_table)

# Save summary table
write.csv(summary_table, file.path(output_dir, "qalb_summary_table.csv"), row.names = FALSE)
cat("\nSummary table saved: qalb_summary_table.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("All results saved to: %s/\n", output_dir))
