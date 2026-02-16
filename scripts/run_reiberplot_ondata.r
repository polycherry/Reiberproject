source("scripts/igg_reiberplot.r")

# Load required library
library(readxl)

# Load actual data from Excel file
data_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/Merged_CSF_Serum_Table_CORRECTED.xlsx"
raw_data <- read_excel(data_path)

# Load age metadata
age_path <- "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/raw_data/metadata_subjectID_age.csv"
age_data <- read.csv(age_path)

# Prepare data frame for Reiber plot
# Convert column names and prepare in the format expected by plot_reibergram
dfa <- data.frame(
  patient_id = raw_data$Sample,
  group = raw_data$Group,
  type = raw_data$Type,
  csf_igg = as.numeric(raw_data$CSF_IgG),
  s_igg = as.numeric(raw_data$Serum_IgG),
  csf_alb = as.numeric(raw_data$CSF_Albumin),
  s_alb = as.numeric(raw_data$Serum_Albumin)
)

# Merge age data (assuming first column is name/ID, second is age)
colnames(age_data) <- c("patient_id", "age")
dfa <- merge(dfa, age_data, by = "patient_id", all.x = TRUE)

# Convert to g/L (assuming data is in mg/L)
dfa$csf_igg = dfa$csf_igg / 1000
dfa$csf_alb = dfa$csf_alb / 1000

# Remove any rows with missing values
dfa <- dfa[complete.cases(dfa[, c("csf_igg", "s_igg", "csf_alb", "s_alb")]), ]

# Plot with different colors for Case and Control
# subset_type can be: "full", "baseline", or "followup"
reiber = plot_reibergram(dfa, color_by_group = TRUE, subset_type = "baseline")
ggsave("results/results_reiber.pdf", plot = reiber)
