# IgG Reibergram Plotting Function

This repository contains an R function `plot_reibergram()` that generates **IgG Reibergrams** for a cohort of patients based on cerebrospinal fluid (CSF) and serum measurements. The Reibergram is a widely used tool in neuroimmunology to visualize **intrathecal IgG synthesis** and **blood-CSF barrier function**.

---

## Table of Contents

- [Requirements](#requirements)  
- [Function: `plot_reibergram`](#function-plot_reibergram)  
- [Input Data](#input-data)  
- [Example Usage](#example-usage)  
- [Output](#output)  
- [Notes](#notes)  
- [References](#references)  

---

## Requirements

- R (≥3.6 recommended)  
- Packages:  

```r
library(ggplot2)
library(reshape2)


Function: plot_reibergram
R

plot_reibergram(patient_df)
Purpose: Plot a Reibergram with physiological curves, intrathecal fraction lines, and patient points labeled by ID.

Returns: A ggplot object.


Input Data
The input patient_df must be a data frame with the following mandatory columns:

Column	Type	Description
patient_id	character or factor	Unique identifier for each patient
age	numeric	Patient age in years (used for QAlb cutoff)
csf_alb	numeric	CSF albumin concentration (same units as serum albumin)
s_alb	numeric	Serum albumin concentration
csf_igg	numeric	CSF IgG concentration (same units as serum IgG)
s_igg	numeric	Serum IgG concentration

Markdown

# IgG Reibergram Plotting Function

This repository contains an R function `plot_reibergram()` that generates **IgG Reibergrams** for a cohort of patients based on cerebrospinal fluid (CSF) and serum measurements. The Reibergram is a widely used tool in neuroimmunology to visualize **intrathecal IgG synthesis** and **blood-CSF barrier function**.

***

## Table of Contents

- [Requirements](#requirements)
- [Function: `plot_reibergram`](#function-plot_reibergram)
- [Input Data](#input-data)
- [Example Usage](#example-usage)
- [Output](#output)
- [Notes](#notes)
- [References](#references)

***

## Requirements

- R ($\ge$3.6 recommended)
- Packages:

```r
library(ggplot2)
library(reshape2)
Function: plot_reibergram
R

plot_reibergram(patient_df)
Purpose: Plot a Reibergram with physiological curves, intrathecal fraction lines, and patient points labeled by ID.

Returns: A ggplot object.

Input Data
The input patient_df must be a data frame with the following mandatory columns:

Column	Type	Description
patient_id	character or factor	Unique identifier for each patient
age	numeric	Patient age in years (used for QAlb cutoff)
csf_alb	numeric	CSF albumin concentration (same units as serum albumin)
s_alb	numeric	Serum albumin concentration
csf_igg	numeric	CSF IgG concentration (same units as serum IgG)
s_igg	numeric	Serum IgG concentration

Export to Sheets
Notes:

Units of CSF and serum proteins must be consistent.

Additional columns are allowed but ignored by the function.

The function computes QAlb = csf_alb/s_alb and QIgG =  csf_igg/s_igg internally

Output
A log-log plot of QIgG vs. QAlb with: Physiological curves (Qmean, Qlim, Qlow)
Intrathecal fraction lines (20%, 40%, 60%, 80%)
Vertical line at age-dependent QAlb cutoff
Patient points with IDs
The function returns a ggplot object, which can be further customized if needed.


Notes
The function is designed to not modify the data except for computing QAlb and QIgG
Negative or zero values are set to NA for plotting on a log scale.
The color palette and line types are predefined for clarity but can be changed by modifying the function.
Use a cohort-level maximum age to compute QAlb cutoff; individual age adjustments are optional but not implemented.

Formulas


\[
Q_{\text{lim}}(\text{Ig}) = a \cdot \sqrt{Q_{\text{Alb}}^2 + b^2} - c
\]

\[
Q_{\text{Alb, upper}} = \left( 4 + \frac{\text{age (years)}}{15} \right) \times 10^{-3}
\]

The Ig quotient for a given intrathecal fraction \(f_\text{IF}\):

\[
Q_\text{Ig} = \frac{Q_\text{lim}}{1 - f_\text{IF}}
\]



References
Reiber H, Peter JB. Cerebrospinal fluid analysis: Disease-related data patterns and evaluation programs. Clin Chem 2001;47:903–916.
Reiber H. Dynamics of brain-derived proteins in cerebrospinal fluid. Clin Chim Acta 1994; 232: 11–26.

Author: Vasishta Polisetty
Date: 2025