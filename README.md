# GeoMx Modular Spatial Transcriptomics Analysis Pipeline

This pipeline provides a robust, modular, and reproducible workflow for multi-level differential expression analysis of GeoMx spatial transcriptomics data. It is designed for clarity, extensibility, and ease of use, with comprehensive validation and diagnostic checks at each stage.

## Overview

The pipeline performs:
- SCE annotation and validation
- Processing and normalization
- Pathway-level contrasts from GSVA output
- Spatial deconvolution and cell type scoring
- Transcript-level contrasts with MAST
- Cell proportion heatmaps and statistical analysis
- Volcano plots for differential expression
- Post hoc analyses:
  - Composite myeloid GSVA barplot (top 30 myeloid samples)
  - Transcript-level volcano plots for top 5 Reactome pathways
  - Top 20 and top 100 transcript volcano plots for MAST
- Summary statistics and CSV outputs

## Prerequisites

- R with required packages: `readr`, `dplyr`, `ggrepel`, `pheatmap`, `tibble`, `ggplot2`, `tidyr`, `httr`
- Output files from single-cell pipeline:
  - GSVA pathway scores
  - Spatial deconvolution cell type proportions
  - MAST differential expression results
- Sufficient disk space for output files
- (Optional for post hoc): Clone the gene set repository:
  ```sh
  git clone https://github.com/andrewrech/single.cell.git
  ```

## Quick Start

1. **Clone this repository**:
   ```sh
   git clone https://github.com/nakialw/geomx-analysis-pipeline.git
   cd geomx-analysis-pipeline
   ```

2. **Review the template**: Open `2.analysis/geomx-analysis-pipeline_template.md`

3. **Modify parameters**: Update the user parameters section for your project

4. **Run the analysis**: Copy the R code from the template and execute in R

## Input Data Requirements

### Required Files
1. **GSVA pathway scores** - CSV file with Rct columns and metadata columns
2. **Spatial deconvolution results** - CSV file with cell type proportions
3. **MAST differential expression results** - CSV file with transcript-level statistics

### File Naming Conventions
- GSVA file should contain "gsva_graphics_hm.csv" in the filename
- Spatial deconvolution file should contain "sdctm_sd_beta_hm.csv" in the filename  
- MAST results file should contain "mast_sub_de.csv" in the filename

### Expected Input File Structure
```
data_directory/
├── *gsva_graphics_hm.csv          # GSVA pathway scores
├── *sdctm_sd_beta_hm.csv          # Spatial deconvolution results
└── *mast_sub_de.csv               # MAST differential expression results
```

## Configuration Parameters

### Input/Output Paths
```r
root_dir <- "~/path/to/geomx/data"
project_name <- "geomx_analysis_name"
```

### File Patterns
```r
gsva_file_pattern <- "gsva_graphics_hm.csv"
sdctm_file_pattern <- "sdctm_sd_beta_hm.csv"
mast_file_pattern <- "mast_sub_de.csv"
```

### Analysis Parameters
```r
subject_ref <- "04"  # Baseline subject for contrasts
tcell_threshold <- 55  # Cut-off for "Tcell_high" classification
safeTME_T <- c("CD4.T.cells","CD8.T.cells", "Treg")
```

### Output Options
```r
save_heatmaps <- TRUE
save_csv_results <- TRUE
create_volcano_plot <- TRUE
```

## Workflow Components

### 1. SCE Annotation & Validation
- Annotates and validates SCE objects, ensuring correct structure and metadata.

### 2. Processing & Normalization
- Loads and processes GSVA, spatial deconvolution, and MAST results.
- Adds T-cell scores and groupings.

### 3. Cell Composition Analysis
- **Heatmaps**: Hierarchical clustering by timepoint, subject, region type, and all annotations.
- **Statistical Analysis**: Volcano-style plot showing differential cell proportions between timepoints.

### 4. Pathway-Level Analysis
- **Post vs Pre**: Compares post-infusion vs pre-infusion timepoints.
- **High vs Low T-cell**: Compares high vs low T-cell infiltration regions.
- **Interaction**: Combined effect of timepoint and T-cell status.
- Results saved as CSV with FDR-adjusted p-values.

### 5. Transcript-Level Analysis (MAST)
- Processes MAST differential expression results.
- **Volcano plots**:
  - Top 20 DE transcripts (red, labeled)
  - Top 100 DE transcripts (red, labeled)
- Results saved as CSV and PNG.

### 6. Post Hoc Analyses
- **Composite Myeloid GSVA Barplot**: Computes a composite myeloid score from macrophages, monocytes, mDCs, pDC, and neutrophils. Plots the median GSVA score for selected myeloid pathways in the top 30 myeloid samples.
- **Transcript-level Volcano Plots for Top Pathways**: For the top 5 Reactome pathways (by GSVA), loads gene sets from a local clone of the gene set repository, filters MAST results for those genes, and plots volcano plots (top 20 most significant transcripts for each pathway).

### 7. Summary Statistics
- Reports percentage of significant pathways for each contrast.
- Provides overview of analysis results.

## Output Structure

```
analysis_results/
├── gbm_geomx_I_cell_heatmap_subject.png
├── gbm_geomx_I_cell_heatmap_region.png
├── gbm_geomx_I_cell_heatmap_timepoint.png
├── gbm_geomx_I_cell_heatmap_all_annotations.png
├── gbm_geomx_I_cell_proportion_statistics.png
├── gbm_geomx_I_gsva_post_vs_pre.csv
├── gbm_geomx_I_gsva_high_vs_low.csv
├── gbm_geomx_I_gsva_postHigh_vs_preLow.csv
├── gbm_geomx_I_mast_processed.csv
├── gbm_geomx_I_volcano_top20.png
├── gbm_geomx_I_volcano_top100.png
├── gbm_geomx_I_posthoc_myeloid_gsva_median_barplot.png
├── gbm_geomx_I_volcano_<pathway>.png (for each top Reactome pathway)
└── gbm_geomx_I_summary_stats.csv
```

## Usage

1. **Review the template**: Open `2.analysis/geomx-analysis-pipeline_template.md`
2. **Modify configuration parameters** in the user parameters section
3. **Ensure input files are available** in the specified directory
4. **(Optional for post hoc)**: Clone the gene set repository to your machine
5. **Copy and run the R code** in your R environment
6. **Check output directory** for results and visualizations

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

### Additional Troubleshooting
- **Missing columns**: Verify input files have required metadata columns
- **Gene set file errors**: For post hoc analyses, ensure the gene set files exist locally and the pathway names match the file names
- **Path errors**: Check that all paths are accessible and writable

## Template Features

The analysis pipeline template (`2.analysis/geomx-analysis-pipeline_template.md`) includes:

- **Complete R script** with no placeholders
- **Clear section breaks** with descriptive headers
- **Prerequisites and input requirements**
- **Detailed explanations** for each analysis component
- **Customization guidance** for different projects
- **Troubleshooting section** for common issues
- **Output file descriptions** and organization

## Notes
- The pipeline automatically finds files using pattern matching
- All analyses include progress reporting and debug information
- Results are saved in both CSV and visualization formats
- The script creates necessary output directories automatically
- Post hoc analyses require local access to gene set files
- The template is designed to be copy-paste ready with minimal modifications 