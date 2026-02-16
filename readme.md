# CSF Reibergram Analysis Suite

This repository contains R scripts for comprehensive analysis of cerebrospinal fluid (CSF) and serum measurements, with a focus on **IgG Reibergram visualization**, **intrathecal IgG synthesis**, and **blood-CSF barrier function** assessment in neuroimmunological studies.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Scripts](#scripts)
  - [IgG Reibergram Plotting](#1-igg-reibergram-plotting)
  - [Case vs Control Comparison](#2-case-vs-control-comparison)
  - [Baseline vs Followup Analysis](#3-baseline-vs-followup-analysis)
  - [QAlb Cutoff Analysis](#4-qalb-cutoff-analysis)
- [Input Data Format](#input-data-format)
- [Usage](#usage)
- [Output](#output)
- [References](#references)

---

## Overview

This analysis suite provides tools for:
- **Reibergram visualization** with case/control group coloring
- **Statistical comparison** of CSF/serum biomarkers between groups
- **Longitudinal analysis** of paired baseline-followup measurements
- **Individualized QAlb assessment** against age-adjusted cutoffs

---

## Project Structure

```
Reiberproject/
├── scripts/
│   ├── igg_reiberplot.r              # Core Reibergram plotting function
│   ├── run_reiberplot_ondata.r       # Load data and generate Reibergram
│   ├── compare_case_control.r        # Case vs Control statistical analysis
│   ├── compare_baseline_followup.r   # Paired baseline-followup comparisons
│   └── analyze_qalb_cutoff.r         # Individual QAlb cutoff analysis
├── data/                              # Data files (ignored by git)
├── results/                           # Output plots and tables (ignored by git)
│   ├── baseline_vs_followup/         # Paired comparison plots
│   └── qalb_analysis/                # QAlb analysis plots
└── readme.md
```

---

## Requirements

- **R** (≥3.6 recommended)
- **Required packages**:

```r
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
```

Install missing packages:
```r
install.packages(c("ggplot2", "readxl", "dplyr", "tidyr", "reshape2"))
```

---

## Scripts

### 1. IgG Reibergram Plotting

**Files**: `igg_reiberplot.r`, `run_reiberplot_ondata.r`

**Purpose**: Generate IgG Reibergrams with physiological curves, intrathecal fraction lines, and patient data points.

**Features**:
- Color-coded by group (Case = reddish, Control = blue)
- Subset selection: Full dataset, Baseline only, or Followup only
- Outlier detection and smart labeling (points above Qlim or QAlb cutoff)
- Age-adjusted QAlb cutoff line
- Jittered labels to prevent overlap

**Function signature**:
```r
plot_reibergram(patient_df, color_by_group = TRUE, subset_type = "baseline")
```

**Parameters**:
- `patient_df`: Data frame with required columns
- `color_by_group`: Boolean, color by Case/Control groups
- `subset_type`: "full", "baseline", or "followup"

---

### 2. Case vs Control Comparison

**File**: `compare_case_control.r`

**Purpose**: Statistical comparison of CSF/serum biomarkers between Case and Control groups at baseline.

**Features**:
- Wilcoxon rank-sum test for group comparisons
- Boxplots with jittered points for each variable
- Calculated indices: IgG Index, Albumin Quotient (QAlb), IgG Quotient (QIgG)
- Composite summary plot highlighting significant differences

**Output**:
- Individual plots for each variable (`results/plot_*.png`)
- Composite significance plot (`results/plot_composite_all_comparisons.png`)

**Variables analyzed**:
- CSF Protein, CSF IgG, Serum IgG, CSF Albumin, Serum Albumin, Age
- IgG Index (QIgG/QAlb), Albumin Quotient, IgG Quotient

---

### 3. Baseline vs Followup Analysis

**File**: `compare_baseline_followup.r`

**Purpose**: Paired comparison of measurements between baseline and followup timepoints for patients with both measurements.

**Features**:
- Identifies patients with both baseline and followup data
- Paired Wilcoxon signed-rank test
- Before-after plots with connecting lines for each patient
- Color-coded by Case/Control group
- Summary plot showing significance of changes

**Output**:
- Individual paired plots (`results/baseline_vs_followup/paired_*.png`)
- Summary significance plot (`results/baseline_vs_followup/paired_summary_all_comparisons.png`)

---

### 4. QAlb Cutoff Analysis

**File**: `analyze_qalb_cutoff.r`

**Purpose**: Assess blood-CSF barrier function relative to individualized age-adjusted cutoffs.

**Features**:
- Calculates individual QAlb cutoff: `(4 + age/15) × 0.001`
- Computes QAlb Ratio (observed/expected)
- Identifies patients with barrier dysfunction (ratio > 1)
- Multiple visualization styles

**Output** (all in `results/qalb_analysis/`):
1. `qalb_ratio_boxplot.png` - Boxplot of QAlb/Cutoff ratios
2. `qalb_vs_age_cutoff.png` - QAlb vs age with cutoff curve
3. `qalb_distance_from_cutoff.png` - Absolute distance from cutoff
4. `qalb_ratio_vs_age.png` - Scatter plot of ratios across age
5. `qalb_ratio_violin.png` - Distribution comparison
6. `qalb_summary_table.csv` - Summary statistics by group

---

## Input Data Format

### Required Data Files

1. **Main data**: `Merged_CSF_Serum_Table_CORRECTED.xlsx`
   - Columns: `Sample`, `Group`, `Type`, `CSF_Protein`, `CSF_IgG`, `Serum_IgG`, `CSF_Albumin`, `Serum_Albumin`

2. **Age metadata**: `metadata_subjectID_age.csv`
   - Columns: Patient ID, Age (years)

### Data Frame Structure

Required columns for Reibergram:
- `patient_id`: Unique identifier
- `group`: "Case" or "Control"
- `type`: "Baseline" or "Followup"
- `age`: Patient age in years
- `csf_alb`: CSF albumin (mg/L)
- `s_alb`: Serum albumin (g/L)
- `csf_igg`: CSF IgG (mg/L)
- `s_igg`: Serum IgG (g/L)

**Note**: Units are automatically converted to g/L where needed (CSF values divided by 1000).

---

## Usage

### Generate Reibergram
```r
source("scripts/run_reiberplot_ondata.r")
```

### Run Case vs Control Analysis
```r
source("scripts/compare_case_control.r")
```

### Run Baseline-Followup Analysis
```r
source("scripts/compare_baseline_followup.r")
```

### Run QAlb Cutoff Analysis
```r
source("scripts/analyze_qalb_cutoff.r")
```

---

## Output

All plots are saved to the `results/` directory:
- High resolution PNG files (300 dpi)
- Publication-ready quality
- Consistent color scheme (Cases: #D55E00, Controls: #56B4E9)

---

## Formulas

### Reibergram Curves

**Upper limit of normal (Qlim)**:
$$Q_{\text{lim}}(\text{IgG}) = 0.93 \cdot \sqrt{Q_{\text{Alb}}^2 + 6 \times 10^{-6}} - 1.7 \times 10^{-3}$$

**Mean physiological (Qmean)**:
$$Q_{\text{mean}} = 0.65 \cdot \sqrt{Q_{\text{Alb}}^2 + 8 \times 10^{-6}} - 1.4 \times 10^{-3}$$

**Lower limit (Qlow)**:
$$Q_{\text{low}} = 0.33 \cdot \sqrt{Q_{\text{Alb}}^2 + 2 \times 10^{-6}} - 0.3 \times 10^{-3}$$

### QAlb Cutoff (Blood-CSF Barrier)

$$Q_{\text{Alb, cutoff}} = \left( 4 + \frac{\text{age (years)}}{15} \right) \times 10^{-3}$$

### Calculated Indices

**IgG Index**:
$$\text{IgG Index} = \frac{Q_{\text{IgG}}}{Q_{\text{Alb}}} = \frac{\text{CSF IgG / Serum IgG}}{\text{CSF Albumin / Serum Albumin}}$$

**Albumin Quotient**:
$$Q_{\text{Alb}} = \frac{\text{CSF Albumin}}{\text{Serum Albumin}}$$

**IgG Quotient**:
$$Q_{\text{IgG}} = \frac{\text{CSF IgG}}{\text{Serum IgG}}$$

---

## References

- Reiber H, Peter JB. Cerebrospinal fluid analysis: Disease-related data patterns and evaluation programs. *Clin Chem* 2001;47:903–916.
- Reiber H. Dynamics of brain-derived proteins in cerebrospinal fluid. *Clin Chim Acta* 1994;232:11–26.

---

## Author

Vasishta Polisetty  
Date: February 2026