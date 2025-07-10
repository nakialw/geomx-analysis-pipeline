# GeoMx Modular Analysis Pipeline Template (R)

This template contains the complete R script for the modular GeoMx analysis pipeline. This script performs pathway-level contrasts from GSVA output, spatial deconvolution T-cell scoring, transcript-level contrasts with MAST, and generates comprehensive visualizations and results.

## Overview

The pipeline performs the following analyses:
- **Pathway-level contrasts** from GSVA output (Post vs Pre infusion, T-cell high vs low, interaction)
- **Spatial deconvolution** T-cell score calculation and high/low classification
- **Transcript-level contrasts** with MAST differential expression analysis
- **Cell proportion heatmaps** with hierarchical clustering and annotations
- **Statistical plots** for cell proportions and transcript expression
- **Post hoc analyses** including myeloid composite scoring and pathway-specific transcript analysis

---

## Prerequisites

Before running this pipeline, ensure you have:
- R installed with required packages (see Libraries section)
- Input data files in the correct format
- Sufficient disk space for output files

---

## Input Data Requirements

### Required Files
1. **GSVA pathway scores** - CSV file with Rct columns and metadata columns
2. **Spatial deconvolution results** - CSV file with cell type proportions
3. **MAST differential expression results** - CSV file with transcript-level statistics

### File Naming Conventions
- GSVA file should contain "gsva_graphics_hm.csv" in the filename
- Spatial deconvolution file should contain "sdctm_sd_beta_hm.csv" in the filename  
- MAST results file should contain "mast_sub_de.csv" in the filename

---

## Complete R Script

### Section 0: User Parameters
*Configure these parameters for your specific project*

```r
###############################################################################
# GeoMx Multi-level Differential Expression Pipeline (Modular Version)
# ---------------------------------------------------------------------------
#  - Pathway-level contrasts from GSVA output
#  - Spatial-decon T-cell score & high/low tagging
#  - Transcript-level contrasts with MAST
#  - Writes tidy CSV results for downstream plotting
###############################################################################

#### 0. USER PARAMETERS ######################################################
# Modify these parameters for your specific project

# Input/Output paths
root_dir <- "~/path/to/geomx/data"
project_name <- "geomx_analysis_name"

# File paths (modify these to match your actual file structure)
gsva_file_pattern <- "gsva_graphics_hm.csv"  # Pattern to find GSVA file
sdctm_file_pattern <- "sdctm_sd_beta_hm.csv"  # Pattern to find spatial deconvolution file
mast_file_pattern <- "mast_sub_de.csv"  # Pattern to find MAST results file

# Analysis parameters
subject_ref <- "04"  # Baseline subject for contrasts
tcell_threshold <- 55  # Cut-off for "Tcell_high" classification

# T-cell types for scoring (modify based on your spatial deconvolution output)
safeTME_T <- c("CD4.T.cells","CD8.T.cells", "Treg")

# Output options
save_heatmaps <- TRUE
save_csv_results <- TRUE
create_volcano_plot <- TRUE
```

### Section 1: Load Required Libraries
*Install these packages if not already installed*

```r
###############################################################################

#### 1. LIBRARIES ############################################################
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggrepel)
  library(pheatmap)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(httr)
})
```

### Section 2: Helper Functions
*These functions handle file finding, statistical analysis, and visualization*

```r
#### 2. HELPER FUNCTIONS #####################################################

# Function to find files by pattern
find_file_by_pattern <- function(directory, pattern) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    stop(paste("No files found matching pattern:", pattern, "in directory:", directory))
  }
  if (length(files) > 1) {
    warning(paste("Multiple files found matching pattern:", pattern, ". Using first one:", files[1]))
  }
  return(files[1])
}

# Function to run GLM contrast
run_glm <- function(expr_mat, meta_data, variable, level_ref, level_comp) {
  # step 1: subset metadata and expression matrix to just the two levels
  idx  <- meta_data[[variable]] %in% c(level_ref, level_comp)
  msub <- meta_data[idx, ]
  esub <- expr_mat[, idx, drop = FALSE]
  
  # step 2: make sure the reference level is properly set
  msub[[variable]] <- relevel(factor(msub[[variable]]), ref = level_ref)
  
  # step 3: for each pathway (row), run a GLM comparing group levels
  res <- do.call(rbind, lapply(seq_len(nrow(esub)), function(i) {
    y  <- esub[i, ]
    df <- data.frame(y = as.numeric(y), group = msub[[variable]])
    fit <- glm(y ~ group, data = df)
    coef <- summary(fit)$coefficients[2, ]  # Coefficients for comparison group
    c(logFC = coef[1], pval = coef[4])
  })) |> as_tibble()
  
  # step 4: name the columns properly
  colnames(res) <- c("logFC", "pval")
  
  # step 5: add pathway names, adjust p-values, and order by significance
  res <- res %>%
    mutate(pathway = rownames(esub), .before = 1,
           adj.P.Val = p.adjust(pval, method = "fdr")) %>%
    arrange(adj.P.Val)
  
  return(res)
}

# Function to create cell proportion heatmaps using pheatmap
create_cell_heatmaps <- function(beta_path, output_dir, project_name) {
  cat("Creating cell proportion heatmaps with hierarchical clustering and annotations...\n")
  
  # Load and prepare data
  beta_df <- readr::read_csv(beta_path, show_col_types = FALSE)
  cat("Loaded data with", nrow(beta_df), "rows and", ncol(beta_df), "columns\n")
  
  # Extract metadata
  cell_meta <- beta_df %>% 
    dplyr::select(rn, subject, region_type, timepoint) %>%
    dplyr::mutate(
      subject = factor(as.character(subject)),
      region_type = factor(as.character(region_type)),
      timepoint = factor(as.character(timepoint))
    ) %>%
    as.data.frame()
  rownames(cell_meta) <- cell_meta$rn
  
  # Extract and clean numeric data
  cell_df <- beta_df %>% 
    dplyr::select(-rn, -subject, -region_type, -timepoint, -replicate) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(as.character(.x)))) %>%
    replace(is.na(.), 0)  # Replace any NA with 0
  
  cat("Cell type matrix dimensions:", nrow(cell_df), "x", ncol(cell_df), "\n")
  
  # Convert to matrix
  cell_matrix <- as.matrix(cell_df)
  rownames(cell_matrix) <- beta_df$rn
  
  # Create annotation (samples = columns)
  annotation_all <- cell_meta[rownames(cell_matrix), c("region_type", "subject", "timepoint")]
  
  # Transpose so cell types are rows, samples are columns
  mat <- t(cell_matrix)
  
  # Color palette
  heat_colors <- colorRampPalette(c("lightgray", "#2166AC", "#08306B"))(100)
  
  # 1. Timepoint annotation
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_all["timepoint", drop=FALSE],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    color = heat_colors,
    main = "Cell-type proportions by Timepoint",
    filename = file.path(output_dir, paste0(project_name, "_cell_heatmap_timepoint.png"))
  )
  
  # 2. Subject annotation
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_all["subject", drop=FALSE],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    color = heat_colors,
    main = "Cell-type proportions by Subject",
    filename = file.path(output_dir, paste0(project_name, "_cell_heatmap_subject.png"))
  )
  
  # 3. Region type annotation
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_all["region_type", drop=FALSE],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    color = heat_colors,
    main = "Cell-type proportions by Region Type",
    filename = file.path(output_dir, paste0(project_name, "_cell_heatmap_region.png"))
  )
  
  # 4. All-in-one annotation
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_all,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    color = heat_colors,
    main = "Cell-type proportions (All Annotations)",
    filename = file.path(output_dir, paste0(project_name, "_cell_heatmap_all_annotations.png"))
  )
  
  cat("All heatmaps saved to:", output_dir, "\n")
  return(list(cell_matrix = cell_matrix, cell_meta = cell_meta))
}
```

### Section 3: Main Analysis Pipeline
*This section contains the core analysis workflow*

```r
#### 3. MAIN ANALYSIS PIPELINE ###############################################

# Create output directory
output_dir <- file.path(root_dir, "analysis_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting GeoMx analysis pipeline...\n")
cat("Project:", project_name, "\n")
cat("Output directory:", output_dir, "\n")
```

### Section 4: Load and Process GSVA Data
*Load GSVA pathway scores and prepare for analysis*

```r
#### 4. LOAD GSVA MATRIX & METADATA ##########################################
cat("\nLoading GSVA matrix and metadata...\n")

# Find GSVA file
gsva_path <- find_file_by_pattern(root_dir, gsva_file_pattern)
cat("Found GSVA file:", gsva_path, "\n")

gsva_mat <- read_csv(gsva_path, show_col_types = FALSE)

keep_cols <- c("rn", "region_type", "replicate", "subject", "timepoint")
rct_cols  <- grep("^Rct", names(gsva_mat), value = TRUE)

gsva_meta <- gsva_mat %>% select(all_of(keep_cols))
expr <- gsva_mat %>% select(all_of(rct_cols))

expr_mat <- t(as.matrix(expr))                # rows = pathways, cols = AOIs
rownames(expr_mat) <- rct_cols
colnames(expr_mat) <- gsva_meta$rn

gsva_meta <- gsva_meta %>%
  mutate(across(c(timepoint, region_type, subject), as.factor)) %>%
  mutate(
    region_type = relevel(region_type, ref = "Tumor"),
    subject     = relevel(subject,     ref = subject_ref)
  )
```

### Section 5: T-Cell Analysis and Cell Proportions
*Process spatial deconvolution data and calculate T-cell scores*

```r
#### 5. ADD T-CELL SCORE FROM SPATIAL-DECON ##################################
cat("\nAdding T-cell scores from spatial deconvolution...\n")

# Find spatial deconvolution file
beta_path <- find_file_by_pattern(root_dir, sdctm_file_pattern)
cat("Found spatial deconvolution file:", beta_path, "\n")

sdctm <- read_csv(beta_path, show_col_types = FALSE)
rownames(sdctm) <- sdctm$rn

# Debug: Check available column names
cat("Available columns in spatial deconvolution data:\n")
print(colnames(sdctm))

# Check if our T-cell columns exist
cat("\nChecking for T-cell columns:\n")
for (col in safeTME_T) {
  exists <- col %in% colnames(sdctm)
  cat("  ", col, ":", if(exists) "FOUND" else "NOT FOUND", "\n")
}

# Find T-cell related columns if our expected ones don't exist
if (!all(safeTME_T %in% colnames(sdctm))) {
  cat("\nT-cell columns not found. Searching for T-cell related columns...\n")
  tcell_cols <- grep("T\\.|Tcell|CD4|CD8|Treg", colnames(sdctm), value = TRUE, ignore.case = TRUE)
  cat("Found T-cell related columns:", paste(tcell_cols, collapse = ", "), "\n")
  
  if (length(tcell_cols) > 0) {
    safeTME_T <- tcell_cols
    cat("Updated safeTME_T to:", paste(safeTME_T, collapse = ", "), "\n")
  } else {
    cat("WARNING: No T-cell columns found. Using all numeric columns except metadata.\n")
    # Use all numeric columns except metadata columns
    metadata_cols <- c("rn", "subject", "region_type", "timepoint", "replicate")
    safeTME_T <- setdiff(colnames(sdctm), metadata_cols)
    cat("Using all available cell types:", paste(safeTME_T, collapse = ", "), "\n")
  }
}

# Calculate T-cell scores
cat("\nCalculating T-cell scores using columns:", paste(safeTME_T, collapse = ", "), "\n")
sdctm <- sdctm %>%
  mutate(Tcell_total = rowSums(across(all_of(safeTME_T))),
         Tcell_group = if_else(Tcell_total > tcell_threshold,
                               "Tcell_high", "Tcell_low"))

# Debug: Check T-cell score distribution
cat("\nT-cell score summary:\n")
print(summary(sdctm$Tcell_total))
cat("T-cell group distribution:\n")
print(table(sdctm$Tcell_group))

# Calculate statistical differences for cell proportions
celltype_cols <- setdiff(colnames(sdctm), c("rn", "subject", "region_type", "timepoint", "replicate", "Tcell_total", "Tcell_group"))

# Calculate statistical differences for cell proportions
cell_stats <- do.call(rbind, lapply(celltype_cols, function(cell_type) {
  # Subset data for this cell type
  cell_data <- sdctm[, c("rn", "timepoint", cell_type)]
  colnames(cell_data)[3] <- "proportion"
  
  # Run t-test between pre and post infusion
  pre_data <- cell_data$proportion[cell_data$timepoint == "Pre-inf"]
  post_data <- cell_data$proportion[cell_data$timepoint == "Post-inf"]
  
  if (length(pre_data) > 0 && length(post_data) > 0) {
    test_result <- t.test(post_data, pre_data, paired = FALSE)
    log2FC <- log2(mean(post_data) / mean(pre_data))
    
    data.frame(
      cell_type = cell_type,
      log2FC = log2FC,
      p_value = test_result$p.value,
      neg_log10_p = -log10(test_result$p.value),
      mean_pre = mean(pre_data),
      mean_post = mean(post_data)
    )
  } else {
    data.frame(
      cell_type = cell_type,
      log2FC = 0,
      p_value = 1,
      neg_log10_p = 0,
      mean_pre = 0,
      mean_post = 0
    )
  }
}))

# Filter to finite values for axis limits
finite_log2FC <- cell_stats$log2FC[is.finite(cell_stats$log2FC)]
finite_neglog10p <- cell_stats$neg_log10_p[is.finite(cell_stats$neg_log10_p)]

# Set default limits if all values are NA
x_min <- if(length(finite_log2FC) > 0) min(finite_log2FC) - 0.1 else -1
x_max <- if(length(finite_log2FC) > 0) max(finite_log2FC) + 0.1 else 1
y_max <- if(length(finite_neglog10p) > 0) max(finite_neglog10p) + 0.5 else 1

cell_stat_plot <- ggplot(cell_stats, aes(x = log2FC, y = neg_log10_p)) +
  geom_point(color = "grey40", size = 3) +
  geom_text_repel(aes(label = cell_type),
                  size = 3, max.overlaps = Inf,
                  box.padding = 0.25, segment.size = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "red") +
  labs(
    title = "Differential Cell Proportions (Post vs Pre Infusion)",
    x = "<-- Pre-inf         Post-inf -->\nfold change (log2)",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.05)))

ggsave(
  file.path(output_dir, paste0(project_name, "_cell_proportion_statistics.png")),
  cell_stat_plot, width = 10, height = 8, dpi = 300
)
cat("Statistical plot of cell proportions saved to:", file.path(output_dir, paste0(project_name, "_cell_proportion_statistics.png")), "\n")

# Merge with GSVA metadata
gsva_meta <- gsva_meta %>%
  left_join(sdctm %>% select(rn, Tcell_total, Tcell_group), by = "rn") %>%
  mutate(Tcell_group = factor(Tcell_group,
                              levels = c("Tcell_low","Tcell_high")))
```

### Section 6: Cell Proportion Heatmaps
*Generate hierarchical clustering heatmaps of cell proportions*

```r
#### 6. CREATE CELL PROPORTION HEATMAPS ######################################
cat("\nCreating cell proportion heatmaps...\n")
cell_results <- create_cell_heatmaps(beta_path, output_dir, project_name)
```

### Section 7: Pathway-Level Analysis
*Perform statistical contrasts on GSVA pathway scores*

```r
#### 7. RUN PATHWAY CONTRASTS ################################################
cat("\nRunning pathway contrasts...\n")

# (a) Post-inf vs Pre-inf
res_post_vs_pre <- run_glm(expr_mat, gsva_meta,
                           "timepoint", "Pre-inf", "Post-inf")

# (b) T-cell high vs low
res_high_vs_low <- run_glm(expr_mat, gsva_meta,
                           "Tcell_group", "Tcell_low", "Tcell_high")

# (c) Combined interaction: Post-inf & T-high vs Pre-inf & T-low
gsva_meta <- gsva_meta %>% mutate(tp_tc = interaction(timepoint, Tcell_group, sep = "_"))
res_postHigh_vs_preLow <- run_glm(expr_mat, gsva_meta,
                                  "tp_tc",
                                  "Pre-inf_Tcell_low",
                                  "Post-inf_Tcell_high")

# Save pathway tables
if (save_csv_results) {
  write_csv(res_post_vs_pre, file.path(output_dir, paste0(project_name, "_gsva_post_vs_pre.csv")))
  write_csv(res_high_vs_low, file.path(output_dir, paste0(project_name, "_gsva_high_vs_low.csv")))
  write_csv(res_postHigh_vs_preLow, file.path(output_dir, paste0(project_name, "_gsva_postHigh_vs_preLow.csv")))
}
```

### Section 8: Transcript-Level Analysis
*Process MAST results and create volcano plots*

```r
#### 8. TRANSCRIPT-LEVEL MAST ANALYSIS #######################################
cat("\nProcessing transcript-level MAST results...\n")

# Find MAST results file
mast_path <- find_file_by_pattern(root_dir, mast_file_pattern)
cat("Found MAST results file:", mast_path, "\n")

mast_raw <- read_csv(mast_path, show_col_types = FALSE)

# Process MAST results
mast_raw$signedFC <- mast_raw$FC  

mast <- mast_raw %>%
  filter(is.finite(signedFC)) %>%
  mutate(log2FC = sign(signedFC) * log2(abs(signedFC) + 1e-6),
         neg_log10_p = -log10(pvalue))

# Create volcano plot
if (create_volcano_plot) {
  ## top-100 by significance
  top100 <- mast %>% arrange(pvalue) %>% slice_head(n = 100)
  
  ## top-20 by significance
  top20 <- mast %>% arrange(pvalue) %>% slice_head(n = 20)
  
  ## ten most-significant on each side for top-20
  top_left_20  <- top20 %>% filter(log2FC < 0) %>% slice_head(n = 10)
  top_right_20 <- top20 %>% filter(log2FC > 0) %>% slice_head(n = 10)
  top_labels_20 <- bind_rows(top_left_20, top_right_20)
  
  # Top 20 volcano plot
  volcano_plot_20 <- ggplot(top20, aes(x = log2FC, y = neg_log10_p)) +
    geom_point(alpha = .7, colour = "red") +
    geom_text_repel(data = top_labels_20, aes(label = name),
                    size = 3, max.overlaps = Inf,
                    box.padding = .25, segment.size = .2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Top 20 DE Transcripts (Post vs Pre Infusion)",
      x = "<- Post-inf         Pre-inf ->",
      y = "-log10(p-value)"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0(project_name, "_volcano_top20.png")), 
         volcano_plot_20, width = 10, height = 8, dpi = 300)

  ## ten most-significant on each side for top-100
  top_left_100  <- top100 %>% filter(log2FC < 0) %>% slice_head(n = 10)
  top_right_100 <- top100 %>% filter(log2FC > 0) %>% slice_head(n = 10)
  top_labels_100 <- bind_rows(top_left_100, top_right_100)
  
  # Top 100 volcano plot
  volcano_plot_100 <- ggplot(top100, aes(x = log2FC, y = neg_log10_p)) +
    geom_point(alpha = .7, colour = "red") +
    geom_text_repel(data = top_labels_100, aes(label = name),
                    size = 3, max.overlaps = Inf,
                    box.padding = .25, segment.size = .2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Top 100 DE Transcripts (Post vs Pre Infusion)",
      x = "<- Post-inf         Pre-inf ->",
      y = "-log10(p-value)"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0(project_name, "_volcano_top100.png")), 
         volcano_plot_100, width = 10, height = 8, dpi = 300)
}

# Save processed MAST results
if (save_csv_results) {
  write_csv(mast, file.path(output_dir, paste0(project_name, "_mast_processed.csv")))
}
```

### Section 9: Summary Statistics
*Generate summary statistics for all analyses*

```r
#### 9. SUMMARY STATISTICS ###################################################
cat("\nGenerating summary statistics...\n")

percent_sig <- function(tbl, alpha = 0.05) {
  mean(tbl$adj.P.Val < alpha) * 100
}

summary_stats <- data.frame(
  Comparison = c("Post vs Pre", "High vs Low T-cell", "Interaction"),
  Percent_Significant = c(
    percent_sig(res_post_vs_pre),
    percent_sig(res_high_vs_low),
    percent_sig(res_postHigh_vs_preLow)
  )
)

cat("\nPathways FDR < 0.05:\n")
cat("  post vs pre  :", round(percent_sig(res_post_vs_pre), 1), "%\n")
cat("  high vs low  :", round(percent_sig(res_high_vs_low), 1), "%\n")
cat("  interaction  :", round(percent_sig(res_postHigh_vs_preLow), 1), "%\n")

# Save summary statistics
if (save_csv_results) {
  write_csv(summary_stats, file.path(output_dir, paste0(project_name, "_summary_stats.csv")))
}
```

### Section 10: Completion Summary
*Display completion message and file summary*

```r
#### 10. COMPLETION ##########################################################
cat("\n" , rep("=", 60), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(rep("=", 60), "\n")
cat("Project:", project_name, "\n")
cat("Output directory:", output_dir, "\n")
cat("Files generated:\n")
if (save_heatmaps) cat("  - Cell proportion heatmaps (PNG)\n")
if (save_csv_results) cat("  - Pathway contrast results (CSV)\n")
if (save_csv_results) cat("  - Processed MAST results (CSV)\n")
if (save_csv_results) cat("  - Summary statistics (CSV)\n")
if (create_volcano_plot) cat("  - Volcano plot (PNG)\n")
cat(rep("=", 60), "\n")
```

### Section 11: Post Hoc Analysis - Myeloid Composite Score
*Additional analysis focusing on myeloid cell populations*

```r
#### 11. POST HOC: Composite Myeloid Score and GSVA Median Plot ##############

cat("\nRunning post hoc myeloid composite and GSVA median plot...\n")

# 1. Compute composite myeloid score
myeloid_cols <- c("macrophages", "monocytes", "mDCs", "pDC", "neutrophils")
sdctm$myeloid_score <- rowSums(sdctm[, myeloid_cols], na.rm = TRUE)

# 2. Get top 30 myeloid samples
top_myeloid_rn <- sdctm %>%
  arrange(desc(myeloid_score)) %>%
  slice_head(n = 30) %>%
  pull(rn)

# 3. Subset gsva_mat to those samples
gsva_subset <- gsva_mat %>% filter(rn %in% top_myeloid_rn)

# 4. Pathways of interest
myeloid_pathways <- c(
  "CS_Mac_CXCL10", "CS_Mono_LYPD2", "CS_Mac_MT1H", "CS_Mac_MARK4", "CS_Mac_SPP1_mouse",
  "TREM2_Mac_Tc", "CS_Mac_MT1H_mouse", "CS_Mac_MARK4_mouse", "IL1B_Mo_Tc", "CS_Mac_F13A1_mouse",
  "CS_Mono_CD300LB", "CS_Mac_APOE_mouse", "CS_Mono_PADI4_mouse", "CS_Mac_CCL18", "CS_Mac_CCL18_mouse",
  "FOLR2_Mac_Tc", "CS_Mac_CXCL10_mouse", "IL4I1_Mac_Tc", "CS_Mono_PADI4", "CS_Mac_F13A1",
  "CS_Mono_CCL3L1_mouse", "CS_Mono_CD300LB_mouse", "CS_Mono_CCL3L1", "CS_Mono_LYPD2_mouse",
  "CS_Mac_SPP1", "CS_Mac_APOE"
)

# 5. Compute median GSVA score for each pathway
median_gsva <- sapply(myeloid_pathways, function(pw) {
  if (pw %in% colnames(gsva_subset)) {
    median(gsva_subset[[pw]], na.rm = TRUE)
  } else {
    NA
  }
})

median_gsva_df <- data.frame(
  pathway = myeloid_pathways,
  median_gsva = median_gsva
) %>% filter(!is.na(median_gsva))

# 6. Plot
posthoc_barplot <- ggplot(median_gsva_df, aes(x = median_gsva, y = reorder(pathway, median_gsva))) +
  geom_col(fill = "#2166AC", width = 0.7) +
  labs(
    title = "Median GSVA Score in Top 30 Myeloid Regions",
    x = "Median GSVA Score (Top 30 Myeloid Regions)",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10)
  )

ggsave(
  file.path(output_dir, paste0(project_name, "_posthoc_myeloid_gsva_median_barplot.png")),
  posthoc_barplot, width = 8, height = 10, dpi = 300
)
cat("Post hoc myeloid GSVA median barplot saved to:", file.path(output_dir, paste0(project_name, "_posthoc_myeloid_gsva_median_barplot.png")), "\n")
```

### Section 12: Post Hoc Analysis - Pathway-Specific Transcript Analysis
*Generate transcript-level volcano plots for top pathways*

```r
#### 12. POST HOC: Transcript-level Volcano Plots for Top Reactome Pathways ####

cat("\nRunning post hoc transcript-level volcano plots for top Reactome pathways...\n")

# 1. Get top 5 pathways
top5_pathways <- head(res_post_vs_pre$pathway, 5)

# 2. Path to the local gene sets directory
geneset_local_dir <- "~/single.cell/inst/extdata/testdata/genesets/src/"

geneset_lists <- list()
for (pw in top5_pathways) {
  file_path <- file.path(geneset_local_dir, paste0(pw, ".txt"))
  cat("Trying local file:", file_path, "\n")
  if (!file.exists(file_path)) {
    cat("  Could not find file for pathway:", pw, "\n")
    next
  }
  gene_list <- readLines(file_path)
  gene_list <- trimws(gene_list)
  gene_list <- gene_list[gene_list != ""]
  assign(paste0(pw, "_transcripts"), gene_list)
  geneset_lists[[pw]] <- gene_list
}

# 3. For each pathway, filter MAST results and plot volcano
for (pw in names(geneset_lists)) {
  gene_list <- geneset_lists[[pw]]
  mast_sub <- mast %>% filter(name %in% gene_list)
  
  if (nrow(mast_sub) == 0) {
    cat("  No MAST results for pathway:", pw, "\n")
    next
  }
  
  # Only keep the top 20 most significant transcripts by p-value
  mast_sub_top20 <- mast_sub %>% arrange(pvalue) %>% slice_head(n = 20)
  
  # Top 10 up/down for labeling (from the top 20)
  top_left  <- mast_sub_top20 %>% filter(log2FC < 0) %>% arrange(pvalue) %>% slice_head(n = 10)
  top_right <- mast_sub_top20 %>% filter(log2FC > 0) %>% arrange(pvalue) %>% slice_head(n = 10)
  top_labels <- bind_rows(top_left, top_right)
  
  volcano_pw <- ggplot(mast_sub_top20, aes(x = log2FC, y = neg_log10_p)) +
    geom_point(alpha = .35, color = "grey40") +
    geom_point(data = mast_sub_top20, color = "red") +
    geom_text_repel(data = top_labels, aes(label = name),
                    size = 3, max.overlaps = Inf,
                    box.padding = .25, segment.size = .2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste0("Transcript-level Volcano Plot: ", pw),
      x = "<- Post-inf         Pre-inf ->",
      y = "-log10(p-value)"
    ) +
    theme_minimal()
  
  ggsave(
    file.path(output_dir, paste0(project_name, "_volcano_", pw, ".png")),
    volcano_pw, width = 10, height = 8, dpi = 300
  )
  cat("  Volcano plot saved for pathway:", pw, "\n")
}
```

---

## Output Files Generated

### Core Analysis Outputs
- **Cell proportion heatmaps** (4 PNG files with different annotations)
- **Statistical cell proportion plot** (volcano-style plot)
- **Pathway contrast results** (3 CSV files for different comparisons)
- **Transcript volcano plots** (top 20 and top 100 MAST results)
- **Summary statistics** (CSV with significance percentages)

### Post Hoc Analysis Outputs
- **Myeloid composite score analysis** (barplot of median GSVA scores)
- **Pathway-specific transcript plots** (volcano plots for top 5 pathways)

---

## Analysis Components

### 1. Pathway-Level Analysis
- GLM contrasts for timepoint (Post vs Pre infusion)
- T-cell group comparisons (High vs Low)
- Interaction analysis (Post-inf & T-high vs Pre-inf & T-low)

### 2. Cell Composition Analysis
- Hierarchical clustering heatmaps with multiple annotations
- Statistical testing of cell proportions between timepoints
- T-cell score calculation and classification

### 3. Transcript-Level Analysis
- MAST results processing and visualization
- Volcano plots for top differentially expressed transcripts
- Pathway-specific transcript analysis

### 4. Post Hoc Analyses
- Myeloid cell population scoring
- Pathway-specific transcript-level investigation

---

## Customization Points

### File Paths and Patterns
- Modify `root_dir` and `project_name` for your project
- Adjust file patterns if your files are named differently
- Update `geneset_local_dir` for pathway-specific analysis

### Analysis Parameters
- Set `subject_ref` to your baseline subject
- Adjust `tcell_threshold` based on your data distribution
- Modify `safeTME_T` to match your cell type column names

### Output Options
- Control which outputs are generated with boolean flags
- Customize plot parameters and file naming

---

## Troubleshooting

### Common Issues
1. **File not found errors**: Check file patterns and directory paths
2. **Missing T-cell columns**: The script will automatically search for alternatives
3. **Memory issues**: Consider processing data in chunks for large datasets
4. **Plot generation errors**: Ensure output directory is writable

### Debug Information
The script provides extensive console output to help identify issues:
- File discovery messages
- Column availability checks
- Data summary statistics
- Progress indicators for each analysis step

This template provides a complete, reproducible analysis pipeline for GeoMx spatial transcriptomics data with comprehensive statistical testing and visualization capabilities.