# Single.cell Analysis Template for GeoMx Data

This template provides a generalized framework for running single.cell analysis on GeoMx spatial transcriptomics data.

## Overview

This pipeline combines:
- Quality Control (QC) analysis
- Batch analysis with GSVA and spatial deconvolution
- Experimental comparison with MAST differential expression (optimized settings)

## Prerequisites

- R with required packages: `andrew.rech`, `single.cell`, `magrittr`, `data.table`, `HDF5Array`
- Annotated SingleCellExperiment (SCE) object
- Docker container with custom R environment (if using containerized analysis)

## Configuration Parameters

### Input/Output Paths
```r
# Modify these paths for your specific project
input_sce_path <- "/path/to/your/annotated_sce.RDS"
output_dir <- "/path/to/output/directory"
project_name <- "your_project_name"
```

### Analysis Parameters
```r
# Assay configuration
assay_idx <- 2  # Which assay to use (typically 2 for exprs)
assay_type <- "wta"  # Assay type
test_type <- "mast"  # Statistical test for differential expression

# Experimental design
timepoints <- c("Pre-inf", "Post-inf")  # Your timepoint levels
subjects <- c("04", "10")  # Your subject IDs
response_factor_order <- c("Pre-inf", "Post-inf")  # Reference level first
```

### Analysis Options
```r
# Toggle which analyses to run
run_qc <- TRUE
run_batch_analysis <- TRUE
run_experimental_comparison <- TRUE

# Output options
save_pdf <- TRUE
save_png <- TRUE
save_rds <- TRUE
save_html <- FALSE
```

### Performance Settings
```r
# Parallelism settings
cores_main <- 1L
cores_assay <- 1L
cores_selections <- 1L
gsva_cores <- 1
gsva_inner_cores <- 4
```

## Expected Input Data Structure

Your SCE object should contain:
- `colData(sce)$timepoint`: Timepoint information
- `colData(sce)$subject`: Subject identifiers
- `colData(sce)$region_type`: Region type information (optional)
- Assay data in the specified assay index

## Expected Output Structure

```
output_directory/
├── qc/
│   └── project_name_*.pdf/png
└── analysis/
    ├── project_name_batch_*.csv/pdf/png
    └── project_name_comparison_*.csv/pdf/png
```

## Key Output Files

### Quality Control
- Metadata quality metrics
- Variance explained plots
- Sample distribution barplots

### Batch Analysis
- GSVA pathway scores
- Spatial deconvolution cell type proportions
- Heatmaps and visualizations

### Experimental Comparison (MAST Differential Expression)
- Gene-level differential expression results
- Statistical significance testing
- Spatial deconvolution cell type proportions
- Heatmaps and visualizations
- Optimized MAST settings for differential expression analysis

## Usage

1. **Modify configuration parameters** at the top of the script
2. **Ensure your SCE object is properly annotated** with required metadata
3. **Run the script** in your R environment
4. **Check output directories** for results and visualizations

## Troubleshooting

### Common Issues

1. **Missing metadata columns**: Ensure your SCE object has the required `timepoint` and `subject` columns in `colData`

2. **Assay index errors**: Verify the correct assay index by checking `names(assays(sce))`

3. **Memory issues**: Reduce parallelism settings or process data in smaller chunks

4. **Path errors**: Ensure all input and output paths are accessible and writable

### Error Handling

The script includes comprehensive error handling with `tryCatch()` blocks. Check the console output for specific error messages and ensure:

- All required libraries are installed
- Input files exist and are readable
- Output directories are writable
- Sufficient memory is available

## Customization

### Adding New Experimental Designs

To add new experimental comparisons, modify the `experimental_design` list:

```r
experimental_design <- list(
  Your_Comparison_Name = list(
    selectorFn = function(sce) {
      # Your selection criteria
      (sce[["your_variable"]] %in% your_levels)
    },
    response = "your_response_variable",
    responseFn = your_response_function,
    fixed = "your_fixed_effects",
    random = "your_random_effects"  # or NULL
  )
)
```

### Modifying Analysis Parameters

Adjust analysis parameters by modifying the `sceAnalyses()` function calls:

#### Batch Analysis (GSVA + Spatial Deconvolution)
```r
sceAnalyses(
  assayIdxs = assay_idx,
  assayType = assay_type,
  l = l,
  tabFilter = "Log2",
  gs = TRUE,           # gene set analysis (GSVA)
  heatmap = TRUE,      # generate heatmaps
  umap = FALSE,        # UMAP dimensionality reduction
  sd = TRUE,           # spatial deconvolution
  rf = FALSE,          # random forest
  # ... other parameters
)
```

#### Experimental Comparison (MAST Differential Expression)
```r
sceAnalyses(
  assayIdxs = assay_idx,
  assayType = assay_type,
  l = l,
  test = "mast",       # MAST statistical test
  tabFilter = "Log2",
  gs = FALSE,          # focus on gene-level analysis
  heatmap = TRUE,      # generate heatmaps
  umap = FALSE,        # UMAP dimensionality reduction
  corr = FALSE,        # disable pairwise correlation
  sd = TRUE,           # spatial deconvolution
  de = TRUE,           # differential expression
  rf = FALSE,          # disable random forest
  # ... other parameters
)
```

## MAST Analysis Optimization

The MAST analysis is optimized for differential expression with these settings:

- **`gs = FALSE`**: Focus on gene-level analysis rather than gene set analysis
- **`corr = FALSE`**: Disable pairwise correlation to focus on differential expression
- **`rf = FALSE`**: Disable random forest to streamline analysis
- **`de = TRUE`**: Enable differential expression analysis
- **`sd = TRUE`**: Include spatial deconvolution for cell type context

## Notes

- This template assumes a two-group comparison (e.g., Pre vs Post treatment)
- Modify the experimental design for more complex comparisons
- The script automatically creates necessary output directories
- All analyses include error handling and progress reporting
- MAST analysis is integrated into the experimental comparison for efficiency 