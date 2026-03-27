library(ggplot2)  # For plotting
library(reshape2) # For melting data frames from wide to long format

# ---------------------------
# Function: plot_reibergram
# ---------------------------
# Purpose: Plot an IgG Reibergram for a cohort of patients
# Input: 
#   patient_df: data frame with columns:
#     - patient_id: unique identifier for each patient (used for labeling points)
#     - age: age of each patient in years
#     - csf_alb: cerebrospinal fluid albumin concentration
#     - s_alb: serum albumin concentration
#     - csf_igg: cerebrospinal fluid IgG concentration
#     - s_igg: serum IgG concentration
#   Notes:
#     - All concentration units should be consistent (e.g., g/L)
#     - patient_df can have multiple patients; the function will compute QAlb/QIgG for each
# Output: ggplot object showing Reibergram
# ---------------------------

plot_reibergram <- function(patient_df, color_by_group = FALSE, subset_type = "full") {
  
  # --- Filter data based on subset type ---
  if (subset_type == "baseline" && "type" %in% colnames(patient_df)) {
    patient_df <- patient_df[patient_df$type == "Baseline", ]
    plot_title <- "IgG Reibergram - Baseline"
  } else if (subset_type == "followup" && "type" %in% colnames(patient_df)) {
    patient_df <- patient_df[patient_df$type == "Followup", ]
    plot_title <- "IgG Reibergram - Followup"
  } else {
    plot_title <- "IgG Reibergram - Full Dataset"
  }
  
  # --- Define Reibergram functions ---
  
  # Age for cohort-based cutoff. Using max age of the cohort as Reibergram cutoffs can be age-dependent.
  age <- max(patient_df$age, na.rm = TRUE)  # ignore NA ages
  # If no age column, use default age of 40
  if (is.na(age) || age == -Inf) {
    age <- 40
  }
  
  cat(sprintf("Using age %.1f years for QAlb cutoff calculation\n", age))
  
  # Create a finely spaced QAlb axis (x-axis of Reibergram)
  qalb <- seq(1e-5, 100e-3, length.out = 1000)  # avoids zero for log scale
  
  # Qlim: upper limit of normal IgG quotient
  # Qmean: mean IgG quotient expected physiologically from population data
  # Qlow: lower limit of normal IgG quotient
  qlim  <- 0.93 * sqrt((qalb ^ 2) + 6e-6) - 1.7e-3
  qmean <- 0.65 * sqrt((qalb ^ 2) + 8e-6) - 1.4e-3
  qlow  <- 0.33 * sqrt((qalb ^ 2) + 2e-6) - 0.3e-3
  
  # Intrathecal fraction lines (dashed lines for various percentages of IgG intrathecal synthesis)
  if20 <- qlim / 0.8  # 20% intrathecal fraction
  if40 <- qlim / 0.6  # 40% intrathecal fraction
  if60 <- qlim / 0.4  # 60% intrathecal fraction
  if80 <- qlim / 0.2  # 80% intrathecal fraction
  
  # QAlb cutoff for blood-CSF barrier dysfunction
  qalb_cutoff <- (4 + age/15) * 0.001  # in same units as qalb (e.g., g/L)
  
  # Combine all curves into one data frame for plotting
  df <- data.frame(qalb, qlim, qmean, qlow, if20, if40, if60, if80)
  
  # Convert wide format to long format suitable for ggplot
  df_long <- melt(df, id.vars = "qalb", variable.name = "curve", value.name = "qval")
  
  # Remove negative values for log scale plotting
  df_long$qval[df_long$qval <= 0] <- NA
  
  # Define line types: dashed for intrathecal fraction lines, solid for physiological curves
  df_long$linetype <- ifelse(grepl("if", df_long$curve), "dashed", "solid")
  
  # Define colors for curves
  curve_colors <- c(
    "qlim"  = "#D55E00",  
    "qmean" = "#0072B2",  
    "qlow"  = "#009E73",  
    "if20"  = "gray50",   
    "if40"  = "gray50",
    "if60"  = "gray50",
    "if80"  = "gray50"
  )
  
  # --- Compute QAlb and QIg for patients ---
  # QAlb = CSF albumin / Serum albumin (reflects blood-CSF barrier)
  patient_df$QAlb <- patient_df$csf_alb / patient_df$s_alb
  # QIgG = CSF IgG / Serum IgG (reflects IgG intrathecal synthesis)
  patient_df$QIgG <- patient_df$csf_igg / patient_df$s_igg
  
  # --- Identify outliers (points above Qlim curve OR QAlb above cutoff) ---
  # For each patient point, calculate expected Qlim at their QAlb
  patient_df$qlim_expected <- 0.93 * sqrt((patient_df$QAlb ^ 2) + 6e-6) - 1.7e-3
  patient_df$is_outlier <- (patient_df$QIgG > patient_df$qlim_expected) | (patient_df$QAlb > qalb_cutoff)
  
  cat(sprintf("QAlb cutoff for blood-CSF barrier: %.3f (×10^-3: %.1f)\n", qalb_cutoff, qalb_cutoff*1e3))
  cat(sprintf("Found %d outlier(s)\n", sum(patient_df$is_outlier)))
  
  # Create a subset for labeling (only outliers)
  outliers_df <- patient_df[patient_df$is_outlier, ]
  
  # Add jitter to label positions to avoid overlap
  if (nrow(outliers_df) > 0) {
    set.seed(42)  # for reproducibility
    # Use multiplicative jitter for log scale - multiply by random factor close to 1
    outliers_df$label_x <- outliers_df$QAlb * (1 + rnorm(nrow(outliers_df), 0, 0.15))
    outliers_df$label_y <- outliers_df$QIgG * (1 + rnorm(nrow(outliers_df), 0.2, 0.15))
    # Ensure labels stay within reasonable bounds
    outliers_df$label_x <- pmax(outliers_df$label_x, 2e-3)
    outliers_df$label_y <- pmax(outliers_df$label_y, 0.5e-3)
  }
  
  # --- Determine colors for groups ---
  if (color_by_group && "group" %in% colnames(patient_df)) {
    point_colors <- c("Case" = "#D55E00", "Control" = "#56B4E9")  # Reddish for cases, blue for controls
    patient_df$plot_color <- patient_df$group
    if (nrow(outliers_df) > 0) {
      outliers_df$plot_color <- outliers_df$group
    }
  } else {
    patient_df$plot_color <- "All"
    if (nrow(outliers_df) > 0) {
      outliers_df$plot_color <- "All"
    }
    point_colors <- c("All" = "black")
  }
  
  # --- Plot ---
  p <- ggplot(df_long, aes(x = qalb, y = qval, color = curve, linetype = linetype)) +
    
    # Draw all Reibergram curves
    geom_line(linewidth = 1) +
    
    # Log10 scale for X-axis (QAlb)
    scale_x_log10(
      breaks = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1),
      labels = function(x) formatC(x*1e3, format="f", digits=0) # convert to per mille for readability
    ) +
    
    # Log10 scale for Y-axis (QIgG)
    scale_y_log10(
      breaks = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1),
      labels = function(x) formatC(x*1e3, format="f", digits=0) # per mille
    ) +
    
    # Add vertical line at QAlb cutoff for blood-CSF barrier dysfunction
    geom_vline(xintercept = qalb_cutoff, linetype = "dashed", color = "black") +
    
    # Assign colors to curves with meaningful labels
    scale_color_manual(values = curve_colors,
                       labels = c("Qlim", "Qmean", "Qlow", "IF 20%", "IF 40%", "IF 60%", "IF 80%")) +
    
    # Use pre-defined linetype
    scale_linetype_identity() +
    
    # Limit visible plot region
    coord_cartesian(xlim = c(2e-3, 100e-3), ylim = c(0.5e-3, 100e-3)) +
    
    # Minimal theme with some customizations
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray85"),
          legend.title = element_blank(),
          legend.position = "right",
          axis.line = element_line(color = "black")) +
    
    # Axis and title labels
    labs(x = expression(Q[Alb]~"(×10"^-3*")"),
         y = expression(Q[Ig]~"(×10"^-3*")"),
         title = plot_title)
  
  # Add patient points with colors based on group (with jitter)
  if (color_by_group && "group" %in% colnames(patient_df)) {
    p <- p + 
      geom_point(data = patient_df, aes(x = QAlb, y = QIgG, color = plot_color),
                 inherit.aes = FALSE, size = 2.5, 
                 position = position_jitter(width = 0.02, height = 0.02, seed = 42)) +
      scale_color_manual(
        values = c(curve_colors, point_colors),
        breaks = names(point_colors),
        labels = names(point_colors),
        name = "Group"
      )
  } else {
    p <- p + 
      geom_point(data = patient_df, aes(x = QAlb, y = QIgG),
                 inherit.aes = FALSE, size = 2, color = "black",
                 position = position_jitter(width = 0.02, height = 0.02, seed = 42))
  }
  
  # Add labels only for outlier points with jitter (color matched to group)
  if (nrow(outliers_df) > 0) {
    if (color_by_group && "group" %in% colnames(outliers_df)) {
      p <- p + 
        geom_text(data = outliers_df, aes(x = label_x, y = label_y, label = patient_id, color = plot_color),
                  inherit.aes = FALSE, size = 3, fontface = "bold")
    } else {
      p <- p + 
        geom_text(data = outliers_df, aes(x = label_x, y = label_y, label = patient_id),
                  inherit.aes = FALSE, size = 3, fontface = "bold", color = "black")
    }
  }
  
  return(p)  # return ggplot object
}
