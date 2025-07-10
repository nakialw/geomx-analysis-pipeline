# Single.cell Processing Pipeline Template (R)

This template contains the complete R script for processing GeoMx spatial transcriptomics data using the single.cell package. This script performs quality control, batch analysis with GSVA and spatial deconvolution, and experimental comparison with MAST differential expression analysis.

## Overview

The pipeline performs the following analyses:
- **Quality Control (QC)** analysis with metadata validation and variance exploration
- **Batch analysis** with GSVA pathway analysis and spatial deconvolution
- **Experimental comparison** with MAST differential expression analysis
- **Data verification** and summary statistics
- **Comprehensive error handling** and progress reporting

## Prerequisites

Before running this pipeline, ensure you have:
- R with required packages: `andrew.rech`, `single.cell`, `magrittr`, `data.table`, `HDF5Array`
- Annotated SingleCellExperiment (SCE) object with proper metadata
- Docker container with custom R environment (if using containerized analysis)
- Sufficient disk space for output files

## Input Data Requirements

### Required SCE Object Metadata
Your SCE object should contain:
- `colData(sce)$timepoint`: Timepoint information (e.g., "Pre-inf", "Post-inf")
- `colData(sce)$subject`: Subject identifiers
- `colData(sce)$region_type`: Region type information (optional)
- Assay data in the specified assay index

### Expected Input File Structure
```
data_directory/
├── annotated_sce.RDS              # Annotated SingleCellExperiment object
└── output_directory/              # Will be created automatically
```

## Complete R Script

### Section 0: User Configuration
*Modify these parameters for your specific project*

```r
# Combined GeoMx Analysis Pipeline using single.cell package
# This script combines single.cell analysis and MAST differential expression analysis
# into a single modular pipeline

# =============================================================================
# USER CONFIGURATION - MODIFY THESE PARAMETERS
# =============================================================================

# Input/Output paths
input_sce_path <- "/path/to/your/annotated_sce.RDS"
output_dir <- "/path/to/output/directory"
project_name <- "your_project_name"

# Analysis parameters
assay_idx <- 2  # Which assay to use (typically 2 for exprs)
assay_type <- "wta"  # Assay type
test_type <- "mast"  # Statistical test for differential expression

# Experimental design parameters
timepoints <- c("Pre-inf", "Post-inf")  # Your timepoint levels
subjects <- c("04", "10")  # Your subject IDs
response_factor_order <- c("Pre-inf", "Post-inf")  # Reference level first

# Analysis options
run_qc <- TRUE
run_batch_analysis <- TRUE
run_experimental_comparison <- TRUE

# Output options
save_pdf <- TRUE
save_png <- TRUE
save_rds <- TRUE
save_html <- FALSE

# Parallelism settings
cores_main <- 1L
cores_assay <- 1L
cores_selections <- 1L
gsva_cores <- 1
gsva_inner_cores <- 4
```

### Section 1: Load Required Libraries and Setup
*Install these packages if not already installed*

```r
# =============================================================================
# LIBRARIES AND SETUP
# =============================================================================

# Load required libraries
library(andrew.rech)
library(single.cell)
library(magrittr)
library(data.table)
library(HDF5Array)

# Setup custom R environment in container
.setup_R()

# Create and set output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\nOutput will be saved to: %s\n", normalizePath(output_dir)))

# Set working directory to output directory
setwd(output_dir)

# Setup output parameters
Sys.setenv(PDF_ANALYSIS_OUTPUT = save_pdf)
Sys.setenv(PNG_ANALYSIS_OUTPUT = save_png)
Sys.setenv(RDS_ANALYSIS_OUTPUT = save_rds)
Sys.setenv(HTML_PLOTLY_ANALYSIS_OUTPUT = save_html)

# Setup graphics device
options(device = function(file = NULL, ...) {
  if(is.null(file)) {
    file <- file.path(output_dir, paste0("temp_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  png(filename = file,
      width = 1200,
      height = 800,
      res = 150,
      type = "cairo",
      bg = "white",
      ...)
})

# Set parallelism
.cores(cores_main)
options(gsva.cores = gsva_cores)
options(gsva.inner.cores = gsva_inner_cores)
```

### Section 2: Data Loading and Verification
*Load your SCE object and verify its structure*

```r
# =============================================================================
# DATA LOADING AND VERIFICATION
# =============================================================================

# Load SCE object
l <- list()
l$sce <- readRDS(input_sce_path)

# Print initial data summary
cat("\nInitial data summary:\n")
cat("Number of samples:", ncol(l$sce), "\n")
cat("Number of genes:", nrow(l$sce), "\n")
cat("Timepoints:", paste(unique(colData(l$sce)$timepoint), collapse=", "), "\n")
cat("Subjects:", paste(unique(colData(l$sce)$subject), collapse=", "), "\n")
cat("Available assays:", paste(names(assays(l$sce)), collapse=", "), "\n")

# Print sample distribution
cat("\nSample distribution by timepoint and subject:\n")
print(table(colData(l$sce)$timepoint, colData(l$sce)$subject, dnn=c("Timepoint", "Subject")))
```

### Section 3: Experimental Design Definition
*Define the statistical framework for your comparisons*

```r
# =============================================================================
# EXPERIMENTAL DESIGN DEFINITION
# =============================================================================

# Define response function
response_fn <- function(x) {
  x <- as.character(x)
  x[x %in% timepoints] <- x[x %in% timepoints]
  factor(x, levels = response_factor_order)
}

# Define experimental design
experimental_design <- list(
  Compare_Pre_Post_Infusion = list(
    selectorFn = function(sce) {
      (sce[["timepoint"]] %in% timepoints) &
      (sce[["subject"]] %in% subjects)
    },
    response = "timepoint",
    responseFn = response_fn,
    responseFactorOrder = response_factor_order,
    fixed = "subject",
    random = NULL,
    responseTab = function(sce) {
      dt <- data.table(
        rn = colnames(sce),
        response = sce[["timepoint"]],
        subject = sce[["subject"]]
      )
      dt[, response := response_fn(response)]
      dt
    }
  )
)
```

### Section 4: Quality Control Analysis
*Run comprehensive quality control checks*

```r
# =============================================================================
# QUALITY CONTROL ANALYSIS
# =============================================================================

if (run_qc) {
  cat("\nRunning Quality Control analysis...\n")
  tryCatch({
    dir.create(file.path(output_dir, "qc"), showWarnings = FALSE, recursive = TRUE)
    
    sceQC(l$sce,
      runQCmetadata = TRUE,
      runQCvarExpl = TRUE,
      runQCbarplots = TRUE,
      fn = file.path(output_dir, "qc", project_name)
    )
    cat("QC analysis completed successfully.\n")
  }, error = function(e) {
    cat("Error in QC step:", conditionMessage(e), "\n")
    print(e)
  })
}
```

### Section 5: Batch Analysis
*Perform GSVA pathway analysis and spatial deconvolution*

```r
# =============================================================================
# BATCH ANALYSIS
# =============================================================================

if (run_batch_analysis) {
  cat("\nRunning batch analysis...\n")
  tryCatch({
    dir.create(file.path(output_dir, "analysis"), showWarnings = FALSE, recursive = TRUE)
    
    selectionsAllRegions() %>%
      sceAnalyses(
        assayIdxs = assay_idx,
        assayType = assay_type,
        l = l,
        tabFilter = "Log2",
        gs = TRUE,
        heatmap = TRUE,
        umap = FALSE,
        sd = TRUE,
        rf = FALSE,
        cores = cores_main,
        assayCores = cores_assay,
        selectionsCores = cores_selections,
        fn = file.path(output_dir, "analysis", paste0(project_name, "_batch"))
      )
    cat("Batch analysis completed successfully.\n")
  }, error = function(e) {
    cat("Error in batch analysis:", conditionMessage(e), "\n")
    print(e)
  })
}
```

### Section 6: Experimental Comparison (MAST)
*Perform differential expression analysis using MAST*

```r
# =============================================================================
# EXPERIMENTAL COMPARISON (INCLUDES MAST)
# =============================================================================

if (run_experimental_comparison) {
  cat("\nRunning experimental comparison analysis with MAST...\n")
  tryCatch({
    experimental_design %>%
      .[1] %>%
      sceAnalyses(
        assayIdxs = assay_idx,
        assayType = assay_type,
        l = l,
        test = test_type,
        tabFilter = "Log2",
        gs = FALSE,  # Focus on gene-level analysis, not gene set analysis
        heatmap = TRUE,
        umap = FALSE,
        corr = FALSE,  # Disable pairwise correlation
        sd = TRUE,
        de = TRUE,
        rf = FALSE,  # Disable random forest
        cores = cores_main,
        assayCores = cores_assay,
        selectionsCores = cores_selections,
        fn = file.path(output_dir, "analysis", paste0(project_name, "_comparison"))
      )
    cat("Experimental comparison with MAST completed successfully.\n")
  }, error = function(e) {
    cat("Error in experimental comparison:", conditionMessage(e), "\n")
    print(e)
  })
}
```

### Section 7: Save Results and Summary
*Save the processed data and display completion summary*

```r
# =============================================================================
# SAVE ANNOTATED SCE OBJECT
# =============================================================================

# Save the annotated SCE object
saveRDS(l$sce, input_sce_path)
cat("\nSaved annotated SCE object to", input_sce_path, "\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n" , rep("=", 60), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(rep("=", 60), "\n")
cat("Output directory:", output_dir, "\n")
cat("Project name:", project_name, "\n")
cat("Analyses run:\n")
cat("  - QC:", if(run_qc) "YES" else "NO", "\n")
cat("  - Batch analysis (GSVA + Spatial Deconvolution):", if(run_batch_analysis) "YES" else "NO", "\n")
cat("  - Experimental comparison (MAST differential expression):", if(run_experimental_comparison) "YES" else "NO", "\n")
cat(rep("=", 60), "\n")
```

---

## Output Files Generated

### Quality Control Outputs
- **QC metadata analysis** - Data quality metrics and validation
- **Variance explained plots** - Understanding data variation
- **Sample distribution barplots** - Visualizing sample distribution

### Batch Analysis Outputs
- **GSVA pathway scores** - Pathway-level analysis results
- **Spatial deconvolution results** - Cell type proportion estimates
- **Heatmaps and visualizations** - Data exploration plots

### Experimental Comparison Outputs
- **MAST differential expression results** - Gene-level statistical testing
- **Differential expression heatmaps** - Visualization of DE results
- **Spatial deconvolution integration** - Cell type context for DE analysis

---

## Analysis Components

### 1. Quality Control
- Metadata validation and quality metrics
- Variance exploration and data distribution analysis
- Sample quality assessment

### 2. Batch Analysis
- GSVA pathway enrichment analysis
- Spatial deconvolution for cell type proportions
- Comprehensive data exploration and visualization

### 3. Experimental Comparison
- MAST differential expression analysis
- Statistical testing between experimental conditions
- Integration with spatial deconvolution results

### 4. Data Management
- Automated output directory creation
- Comprehensive error handling
- Progress reporting and completion summaries

---

## Customization Points

### File Paths and Project Settings
- Modify `input_sce_path` to point to your annotated SCE object
- Update `output_dir` and `project_name` for your project
- Adjust `assay_idx` based on your data structure

### Experimental Design
- Update `timepoints` and `subjects` to match your study design
- Modify `response_factor_order` for proper statistical comparison
- Add additional experimental designs as needed

### Analysis Options
- Toggle individual analyses with boolean flags
- Control output formats (PDF, PNG, RDS, HTML)
- Adjust parallelism settings for your system

### Performance Settings
- Modify core settings based on your system capabilities
- Adjust GSVA parallelism for optimal performance
- Control memory usage through core allocation

---

## Troubleshooting

### Common Issues
1. **Missing metadata columns**: Ensure your SCE object has required `timepoint` and `subject` columns
2. **Assay index errors**: Verify the correct assay index by checking `names(assays(sce))`
3. **Memory issues**: Reduce parallelism settings or process data in smaller chunks
4. **Path errors**: Ensure all input and output paths are accessible and writable

### Debug Information
The script provides extensive console output to help identify issues:
- Data summary and verification messages
- Progress indicators for each analysis step
- Error messages with detailed information
- Completion summaries with analysis status

### Error Handling
- Comprehensive `tryCatch()` blocks prevent script crashes
- Detailed error reporting for troubleshooting
- Graceful handling of missing data or configuration issues

---

## MAST Analysis Optimization

The MAST analysis is optimized for differential expression with these settings:

- **`gs = FALSE`**: Focus on gene-level analysis rather than gene set analysis
- **`corr = FALSE`**: Disable pairwise correlation to focus on differential expression
- **`rf = FALSE`**: Disable random forest to streamline analysis
- **`de = TRUE`**: Enable differential expression analysis
- **`sd = TRUE`**: Include spatial deconvolution for cell type context

---

## Notes
- This template assumes a two-group comparison (e.g., Pre vs Post treatment)
- Modify the experimental design for more complex comparisons
- The script automatically creates necessary output directories
- All analyses include error handling and progress reporting
- MAST analysis is integrated into the experimental comparison for efficiency
- The template is designed to be copy-paste ready with minimal modifications

This template provides a complete, reproducible processing pipeline for GeoMx spatial transcriptomics data with comprehensive quality control, batch analysis, and differential expression capabilities. 