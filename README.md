# GeoMx Spatial Transcriptomics Analysis Pipeline

A comprehensive, modular pipeline for analyzing NanoString GeoMx spatial transcriptomics data, from raw data processing to differential expression analysis.

## Overview

This pipeline provides a complete workflow for GeoMx data analysis, organized into logical stages:

- **00.docker-instance/**: Docker setup and file transfer protocols
- **0.sce-annotation/**: SCE object annotation and metadata management
- **1.processing/**: Single.cell analysis and MAST differential expression
- **2.analysis/**: Multi-level differential expression and visualization

## Pipeline Structure

```
gbm-st-code/
├── 00.docker-instance/           # Docker setup and file transfer
│   ├── docker_template.md        # Master Docker setup guide
│   ├── new-instance-setup.md     # New instance initialization
│   └── abi_base.sh              # Base shell script
├── 0.sce-annotation/            # SCE object annotation
│   ├── annotation_template.md    # Generalized annotation guide
│   └── geomx-annotation.md       # Master annotation script
├── 1.processing/                # Data processing and analysis
│   ├── single.cell_combined.r    # Master processing script (includes MAST)
│   └── single.cell_template.md   # Generalized processing guide
├── 2.analysis/                  # Multi-level analysis and visualization
│   ├── geomx-analysis-pipeline_modular.R  # Master analysis pipeline
│   └── geomx-analysis-pipeline_template.md # Generalized analysis guide
└── README.md                    # This file
```

## Quick Start

### 1. Setup Environment
```bash
# Navigate to project directory
cd gbm-st-code

# Review docker setup (if needed)
cat 00.docker-instance/docker_template.md
```

### 2. Annotate SCE Objects (if needed)
```r
# Review annotation template
# Modify 0.sce-annotation/annotation_template.md for your project
# Run annotation script
source("0.sce-annotation/geomx-annotation.md")
```

### 3. Run Processing Pipeline
```r
# Configure parameters in single.cell_combined.r
# Run combined analysis
source("1.processing/single.cell_combined.r")
```

### 4. Run Analysis Pipeline
```r
# Configure parameters in geomx-analysis-pipeline_modular.R
# Run modular analysis
source("2.analysis/geomx-analysis-pipeline_modular.R")
```

## Detailed Usage

### Stage 0: Docker Setup (Optional)

If you need to set up a new Docker environment:

1. **Review the template**: `00.docker-instance/docker_template.md`
2. **Configure paths**: Modify docker configuration for your system
3. **Set up container**: Follow the setup instructions
4. **Transfer files**: Use the file transfer protocols

### Stage 1: SCE Annotation (Optional)

If you need to annotate your SCE objects:

1. **Review the template**: `0.sce-annotation/annotation_template.md`
2. **Configure experimental design**: Define your experimental groups
3. **Run annotation**: Execute the annotation script
4. **Validate results**: Check the Excel output for verification

### Stage 2: Data Processing

The main processing pipeline combines single.cell analysis and MAST differential expression:

1. **Configure parameters**: Modify the user configuration section in `single.cell_combined.r`
2. **Set input/output paths**: Point to your SCE object and desired output directory
3. **Define experimental design**: Configure timepoints, subjects, and analysis parameters
4. **Run analysis**: Execute the script to generate:
   - Quality control reports
   - GSVA pathway scores
   - Spatial deconvolution results
   - MAST differential expression results

### Stage 3: Multi-level Analysis

The analysis pipeline performs comprehensive downstream analysis:

1. **Configure parameters**: Modify the user configuration section in `geomx-analysis-pipeline_modular.R`
2. **Set file patterns**: Configure patterns to find your analysis outputs
3. **Run analysis**: Execute the script to generate:
   - Cell proportion heatmaps
   - Pathway-level contrasts
   - Transcript-level analysis
   - Summary statistics

## Configuration

### Key Parameters

#### Processing Pipeline (`single.cell_combined.r`)
```r
# Input/Output paths
input_sce_path <- "/path/to/your/annotated_sce.RDS"
output_dir <- "/path/to/output/directory"
project_name <- "your_project_name"

# Analysis parameters
assay_idx <- 2  # Which assay to use
assay_type <- "wta"  # Assay type
test_type <- "mast"  # Statistical test

# Experimental design
timepoints <- c("Pre-inf", "Post-inf")
subjects <- c("04", "10")
```

#### Analysis Pipeline (`geomx-analysis-pipeline_modular.R`)
```r
# Input/Output paths
root_dir <- "/path/to/your/data/directory"
project_name <- "your_project_name"

# File patterns
gsva_file_pattern <- "gsva_graphics_hm.csv"
sdctm_file_pattern <- "sdctm_sd_beta_hm.csv"
mast_file_pattern <- "mast_sub_de.csv"

# Analysis parameters
subject_ref <- "04"  # Baseline subject
tcell_threshold <- 55  # T-cell classification threshold
```

## Output Structure

### Processing Pipeline Outputs
```
output_directory/
├── qc/
│   └── project_name_*.pdf/png
├── analysis/
│   ├── project_name_batch_*.csv/pdf/png
│   └── project_name_comparison_*.csv/pdf/png
└── mast_analysis/
    └── project_name_mast_analysis_*.csv/pdf/png
```

### Analysis Pipeline Outputs
```
analysis_results/
├── project_name_cell_heatmap_*.png
├── project_name_gsva_*.csv
├── project_name_mast_processed.csv
├── project_name_volcano_plot.png
└── project_name_summary_stats.csv
```

## Key Features

### Modularity
- Each stage can be run independently
- Configurable parameters for different projects
- Generalized templates for reusability

### Error Handling
- Comprehensive error checking and reporting
- Graceful failure handling
- Progress reporting throughout

### Flexibility
- Support for different experimental designs
- Configurable analysis parameters
- Multiple output formats

### Reproducibility
- Version-controlled scripts
- Documented parameter configurations
- Generalized templates for future use

## Troubleshooting

### Common Issues

1. **File not found errors**
   - Check file paths and patterns
   - Verify file naming conventions
   - Ensure files are accessible

2. **Memory issues**
   - Reduce parallelism settings
   - Process data in smaller chunks
   - Increase container memory if using Docker

3. **Package dependency issues**
   - Install required R packages
   - Check Bioconductor package versions
   - Verify R environment setup

### Getting Help

1. **Check templates**: Review the `*_template.md` files for guidance
2. **Validate inputs**: Ensure your data matches expected formats
3. **Check logs**: Review console output for specific error messages
4. **Test with subset**: Try running with a smaller data subset first

## Customization

### Adding New Analyses

1. **Modify experimental design**: Update the experimental design configuration
2. **Add new contrasts**: Extend the analysis pipeline with new comparisons
3. **Customize visualizations**: Modify plotting parameters and color schemes
4. **Add new outputs**: Extend the output generation functions

### Adapting for Different Data Types

1. **Update file patterns**: Modify file discovery patterns for your data
2. **Adjust parameters**: Configure analysis parameters for your experimental design
3. **Modify cell types**: Update T-cell type definitions for your spatial deconvolution
4. **Customize thresholds**: Adjust classification thresholds based on your data

## Best Practices

1. **Version control**: Keep scripts under version control
2. **Documentation**: Document any custom modifications
3. **Testing**: Test with small datasets before full analysis
4. **Backup**: Keep backup copies of important data and results
5. **Validation**: Always validate results and check for expected patterns

## Dependencies

### R Packages
- `andrew.rech`
- `single.cell`
- `magrittr`
- `data.table`
- `HDF5Array`
- `readr`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `pheatmap`
- `tibble`
- `digest`
- `openxlsx`

### System Requirements
- R 4.0+
- Sufficient memory for data processing
- Docker (optional, for containerized analysis)

## Contributing

When modifying the pipeline:

1. **Update templates**: Keep generalized templates current
2. **Document changes**: Add comments and update documentation
3. **Test thoroughly**: Ensure modifications work with different datasets
4. **Maintain compatibility**: Preserve backward compatibility where possible

## License

This pipeline is provided for research use. Please cite appropriate references for the underlying methods and tools used.

## Support

For questions or issues:

1. Check the troubleshooting section
2. Review the template files for guidance
3. Validate your data and configuration
4. Test with smaller datasets first

---

**Note**: This pipeline is designed for GeoMx spatial transcriptomics data. Adaptations may be needed for other spatial transcriptomics platforms or experimental designs. 