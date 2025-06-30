# GeoMx Analysis Pipeline Template

This template provides a generalized framework for running multi-level differential expression analysis on GeoMx spatial transcriptomics data.

## Overview

This pipeline performs:
- Pathway-level contrasts from GSVA output
- Spatial deconvolution T-cell scoring and classification
- Transcript-level contrasts with MAST
- Cell proportion heatmaps
- Volcano plots for differential expression
- Summary statistics and CSV outputs

## Prerequisites

- R with required packages: `readr`, `dplyr`, `ggrepel`, `pheatmap`, `tibble`, `ggplot2`
- Output files from single.cell pipeline:
  - GSVA pathway scores
  - Spatial deconvolution cell type proportions
  - MAST differential expression results

## Configuration Parameters

### Input/Output Paths
```r
# Modify these paths for your specific project
root_dir <- "~/path/to/your/data/directory"
project_name <- "your_project_name"
```

### File Patterns
```r
# File patterns to find your analysis outputs
gsva_file_pattern <- "gsva_graphics_hm.csv"
sdctm_file_pattern <- "sdctm_sd_beta_hm.csv"
mast_file_pattern <- "mast_sub_de.csv"
```

### Analysis Parameters
```r
# Analysis configuration
subject_ref <- "04"  # Baseline subject for contrasts
tcell_threshold <- 55  # Cut-off for "Tcell_high" classification

# T-cell types for scoring (modify based on your spatial deconvolution output)
safeTME_T <- c("T.CD4.naive", "T.CD4.memory",
               "T.CD8.naive", "T.CD8.memory", "Treg")
```

### Output Options
```r
# Toggle which outputs to generate
save_heatmaps <- TRUE
save_csv_results <- TRUE
create_volcano_plot <- TRUE
```

## Expected Input File Structure

Your data directory should contain files matching these patterns:

```
data_directory/
├── *gsva_graphics_hm.csv          # GSVA pathway scores
├── *sdctm_sd_beta_hm.csv          # Spatial deconvolution results
└── *mast_sub_de.csv               # MAST differential expression results
```

### GSVA File Format
Expected columns:
- `rn`: Region/sample identifiers
- `region_type`: Region type (e.g., "Tumor", "Stroma")
- `replicate`: Replicate information
- `subject`: Subject identifiers
- `timepoint`: Timepoint information
- `Rct*`: Reactome pathway scores (columns starting with "Rct")

### Spatial Deconvolution File Format
Expected columns:
- `rn`: Region/sample identifiers
- `subject`, `region_type`, `timepoint`: Metadata columns
- Cell type proportion columns (e.g., `T.CD4.naive`, `T.CD8.memory`, etc.)

### MAST Results File Format
Expected columns:
- `name`: Gene/transcript names
- `FC`: Fold change values
- `pvalue`: P-values

## Expected Output Structure

```
analysis_results/
├── project_name_cell_heatmap_subject.png
├── project_name_cell_heatmap_region.png
├── project_name_cell_heatmap_timepoint.png
├── project_name_gsva_post_vs_pre.csv
├── project_name_gsva_high_vs_low.csv
├── project_name_gsva_postHigh_vs_preLow.csv
├── project_name_mast_processed.csv
├── project_name_volcano_plot.png
└── project_name_summary_stats.csv
```

## Analysis Components

### 1. Cell Proportion Heatmaps
- Creates heatmaps of cell type proportions
- Annotated by subject, region type, and timepoint
- Uses hierarchical clustering for visualization

### 2. Pathway-Level Contrasts
- **Post vs Pre**: Compares post-infusion vs pre-infusion timepoints
- **High vs Low T-cell**: Compares high vs low T-cell infiltration regions
- **Interaction**: Combined effect of timepoint and T-cell status

### 3. Transcript-Level Analysis
- Processes MAST differential expression results
- Creates volcano plots highlighting top differentially expressed genes
- Calculates log2 fold changes and significance metrics

### 4. Summary Statistics
- Reports percentage of significant pathways for each contrast
- Provides overview of analysis results

## Usage

1. **Modify configuration parameters** at the top of the script
2. **Ensure input files are available** in the specified directory
3. **Run the script** in your R environment
4. **Check output directory** for results and visualizations

## Customization

### Modifying T-cell Types
Adjust the `safeTME_T` vector to match your spatial deconvolution output:

```r
# Example for different cell type sets
safeTME_T <- c("T.CD4.naive", "T.CD4.memory", "T.CD8.naive", "T.CD8.memory", "Treg")
# or
safeTME_T <- c("CD4_T", "CD8_T", "Treg", "NK")  # Different naming convention
```

### Adding New Contrasts
To add new pathway contrasts, modify the analysis section:

```r
# Example: Adding a new contrast
res_new_contrast <- run_glm(expr_mat, gsva_meta,
                           "your_variable", "reference_level", "comparison_level")

# Save the new results
write_csv(res_new_contrast, file.path(output_dir, paste0(project_name, "_new_contrast.csv")))
```

### Modifying Heatmap Parameters
Adjust heatmap creation in the `create_cell_heatmaps` function:

```r
# Example: Changing color scheme
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Example: Changing clustering method
pheatmap(
  cell_matrix,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  # ... other parameters
)
```

## Troubleshooting

### Common Issues

1. **File not found errors**: Ensure file patterns match your actual file names
2. **Missing columns**: Verify input files have required metadata columns
3. **Memory issues**: Process data in smaller chunks if needed
4. **Path errors**: Check that all paths are accessible and writable

### Error Handling

The script includes error handling for:
- File discovery and loading
- Data processing and analysis
- Output generation and saving

Check console output for specific error messages and ensure:
- All required libraries are installed
- Input files exist and are readable
- Output directories are writable
- Data formats match expected structure

## Notes

- The pipeline automatically finds files using pattern matching
- All analyses include progress reporting
- Results are saved in both CSV and visualization formats
- The script creates necessary output directories automatically
- T-cell threshold should be determined based on your data distribution 