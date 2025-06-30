# SCE Annotation Template for GeoMx Data

This template provides a generalized framework for annotating SingleCellExperiment (SCE) objects with experimental metadata, including comprehensive validation and diagnostic checks.

## Overview

This pipeline performs:
- Loading raw SCE objects
- Adding experimental metadata (timepoints, subjects, etc.)
- Comprehensive validation and sanity checks
- Diagnostic analysis of SCE object structure
- Saving annotated SCE objects

## Prerequisites

- R with required packages: `digest`, `openxlsx`, `data.table`, `SingleCellExperiment`, `andrew.rech`, `single.cell`, `magrittr`, `HDF5Array`
- Raw SCE objects from GeoMx processing
- Experimental metadata information

## Configuration Parameters

### Input/Output Paths
```r
# Modify these paths for your specific project
cosmx_path <- "/path/to/your/cosmx_sce.RDS"
geomx_raw_path <- "/path/to/your/raw_geomx_sce.RDS"
geomx_annotated_path <- "/path/to/your/annotated_geomx_sce.RDS"
output_path <- "/path/to/output/annotated_sce.RDS"
excel_output <- "metadata_sanity_check.xlsx"
```

### Experimental Design
```r
# Define your experimental groups and metadata
experimental_groups <- list(
  group1 = list(
    accession = "SH-23-0028491 A-2",
    timepoint = "Post-inf",
    subject = "04",
    days_since_infusion = -50,
    reinfusion = TRUE
  ),
  group2 = list(
    accession = c("SH-23-0018172 B-2", "SH-23-0018172 A-3"),
    timepoint = "Pre-inf",
    subject = "04",
    days_since_infusion = -30,
    reinfusion = TRUE
  )
  # Add more groups as needed
)
```

### Sentinel Columns
```r
# Columns to use for validation (should not change during annotation)
sentinel_cols <- c("roi", "GenesDetected")
```

## Expected Input Data Structure

Your SCE objects should contain:
- `colData(sce)$accession`: Accession identifiers
- `colData(sce)$slide_name`: Slide identifiers
- `colData(sce)$roi`: Region of interest identifiers
- `colData(sce)$GenesDetected`: Gene detection counts

## Annotation Process

### 1. Load Raw Data
```r
# Load the raw SCE objects
cosmx <- readRDS(cosmx_path)
geomx_raw <- readRDS(geomx_raw_path)

# Check available metadata columns
colnames(colData(cosmx))
colnames(colData(geomx_raw))
```

### 2. Validate Data Structure
```r
# Show unique tissue identifiers and their associated metadata
table(colData(cosmx)$Run_Tissue_name, 
      colData(cosmx)$subject, 
      colData(cosmx)$timepoint)

# Show unique slide identifiers and their associated metadata
table(colData(geomx_raw)$slide_name, 
      colData(geomx_raw)$subject, 
      colData(geomx_raw)$timepoint)
```

### 3. Prepare for Annotation
```r
# Convert to data.table for easier manipulation
cd <- as.data.table(colData(geomx_raw))

# Create sentinel snapshot for validation
cd[, row_id := .I]
before <- cd[, ..sentinel_cols]
hash_before <- digest(unname(as.matrix(before)))
```

### 4. Apply Annotations
```r
# Apply experimental group annotations
for (group_name in names(experimental_groups)) {
  group <- experimental_groups[[group_name]]
  
  if (length(group$accession) == 1) {
    # Single accession
    cd[accession == group$accession,
       `:=`(timepoint = group$timepoint,
            subject = group$subject,
            days_since_infusion = group$days_since_infusion,
            reinfusion = group$reinfusion)]
  } else {
    # Multiple accessions
    cd[accession %chin% group$accession,
       `:=`(timepoint = group$timepoint,
            subject = group$subject,
            days_since_infusion = group$days_since_infusion,
            reinfusion = group$reinfusion)]
  }
}
```

### 5. Validation Checks
```r
# Verify row order hasn't changed
stopifnot(identical(cd$row_id, seq_len(nrow(cd))))

# Verify sentinel values haven't changed
after <- cd[, ..sentinel_cols]
hash_after <- digest(unname(as.matrix(after)))
stopifnot(identical(hash_before, hash_after))

# Clean up helper column
cd[, row_id := NULL]
```

### 6. Save Results
```r
# Update SCE object
colData(geomx_raw) <- DataFrame(cd)

# Create Excel summary
wb <- createWorkbook()
addWorksheet(wb, "geomx_meta")
writeData(wb, "geomx_meta",
          cd[, .(slide_name, accession, subject,
                 timepoint, roi, GenesDetected)][order(slide_name, roi)])
saveWorkbook(wb, excel_output, overwrite = TRUE)

# Save annotated SCE object
saveRDS(geomx_raw, output_path)
```

## Comprehensive Validation and Diagnostic Checks

### 7. Sanity Check Validation
```r
# Load both raw and annotated objects for comparison
sce_raw <- readRDS(geomx_raw_path)
sce_annot <- readRDS(output_path)

# Fix column names if needed (remove .dcc suffix)
ids_raw <- sub("\\.dcc$", "", colnames(sce_raw))
ids_annot <- colnames(sce_annot)

# Replace names in the raw object
colnames(sce_raw) <- ids_raw

# Check identity and order
cat("Identical order now? ", identical(ids_raw, ids_annot), "\n")
cat("Same set now?        ", setequal(ids_raw, ids_annot), "\n")

if (!identical(ids_raw, ids_annot)) {
  cat("First 10 mismatches (raw vs annotated):\n")
  print(head(data.table(raw = ids_raw, annot = ids_annot)[raw != annot], 10))
}

# Sentinel hash check (ROI + GenesDetected)
sentinel <- c("roi", "GenesDetected")
hash_raw <- digest(as.matrix(as.data.table(colData(sce_raw)[, sentinel])))
hash_annot <- digest(as.matrix(as.data.table(colData(sce_annot)[, sentinel])))

cat("Metadata sentinels identical? ", identical(hash_raw, hash_annot), "\n")

# Assay-level equality check
if (setequal(ids_raw, ids_annot)) {
  # Align columns if objects were shuffled
  sce_raw <- sce_raw[, ids_annot]  # now same order
  assay_ok <- identical(assay(sce_raw, 2), assay(sce_annot, 2))
  cat("Assay matrix (idx 2) identical? ", assay_ok, "\n")
  
  if (!assay_ok) {
    diff <- which(assay(sce_raw, 2) != assay(sce_annot, 2), arr.ind = TRUE)
    cat("Non-matching values: ", nrow(diff), "\n")
  }
}

# Report new columns added
new_cols <- setdiff(colnames(colData(sce_annot)), colnames(colData(sce_raw)))
cat("New metadata columns in annotated SCE:\n")
print(new_cols)
```

### 8. Diagnostic Analysis
```r
# Load the annotated SCE object for diagnostics
l <- list()
l$sce <- readRDS(output_path)

# Check basic structure
cat("\n=== SCE Object Structure ===\n")
cat("Class:", class(l$sce), "\n")
cat("Dimensions:", dim(l$sce), "\n")
cat("Number of assays:", length(assays(l$sce)), "\n")
cat("Assay names:", names(assays(l$sce)), "\n")

# Check column data
cat("\n=== Column Data (colData) ===\n")
cat("Number of columns in colData:", ncol(colData(l$sce)), "\n")
cat("Column names in colData:\n")
print(colnames(colData(l$sce)))

# Check specific columns we need
cat("\n=== Checking Required Columns ===\n")
cat("Has 'timepoint' column:", "timepoint" %in% colnames(colData(l$sce)), "\n")
cat("Has 'subject' column:", "subject" %in% colnames(colData(l$sce)), "\n")
cat("Has 'SampleID' column:", "SampleID" %in% colnames(colData(l$sce)), "\n")

# Check column names of SCE object
cat("\n=== SCE Column Names ===\n")
cat("Has column names:", !is.null(colnames(l$sce)), "\n")
if(!is.null(colnames(l$sce))) {
  cat("Number of column names:", length(colnames(l$sce)), "\n")
  cat("First few column names:", paste(head(colnames(l$sce)), collapse=", "), "\n")
} else {
  cat("Column names are NULL\n")
}

# Check timepoint values
if("timepoint" %in% colnames(colData(l$sce))) {
  cat("\n=== Timepoint Values ===\n")
  cat("Unique timepoints:", paste(unique(colData(l$sce)$timepoint), collapse=", "), "\n")
  cat("Timepoint class:", class(colData(l$sce)$timepoint), "\n")
}

# Check subject values
if("subject" %in% colnames(colData(l$sce))) {
  cat("\n=== Subject Values ===\n")
  cat("Unique subjects:", paste(unique(colData(l$sce)$subject), collapse=", "), "\n")
  cat("Subject class:", class(colData(l$sce)$subject), "\n")
}

# Test the selector function
cat("\n=== Testing Selector Function ===\n")
tryCatch({
  selector_result <- (colData(l$sce)$timepoint %in% c("Pre-inf", "Post-inf")) &
                    (colData(l$sce)$subject %in% c("04", "10"))
  cat("Selector function result length:", length(selector_result), "\n")
  cat("Number of TRUE values:", sum(selector_result), "\n")
}, error = function(e) {
  cat("Error in selector function:", conditionMessage(e), "\n")
})
```

### 9. Critical SCE Object Structure Checks
```r
# =============================================================================
# CRITICAL SCE OBJECT STRUCTURE CHECKS
# =============================================================================

cat("\n" , rep("=", 60), "\n")
cat("CRITICAL SCE OBJECT STRUCTURE CHECKS\n")
cat(rep("=", 60), "\n")

# Check 1: rownames(colData) must equal colnames(SCE)
cat("\n=== Check 1: rownames(colData) == colnames(SCE) ===\n")
if (!is.null(rownames(colData(l$sce))) && !is.null(colnames(l$sce))) {
  rownames_colData <- rownames(colData(l$sce))
  colnames_sce <- colnames(l$sce)
  
  cat("rownames(colData) length:", length(rownames_colData), "\n")
  cat("colnames(SCE) length:", length(colnames_sce), "\n")
  cat("Identical order:", identical(rownames_colData, colnames_sce), "\n")
  cat("Same set:", setequal(rownames_colData, colnames_sce), "\n")
  
  if (!identical(rownames_colData, colnames_sce)) {
    cat("WARNING: rownames(colData) != colnames(SCE)\n")
    cat("First 5 mismatches:\n")
    mismatches <- data.frame(
      rownames_colData = head(rownames_colData, 5),
      colnames_sce = head(colnames_sce, 5)
    )
    print(mismatches)
  } else {
    cat("✓ rownames(colData) == colnames(SCE) - PASSED\n")
  }
} else {
  cat("ERROR: Either rownames(colData) or colnames(SCE) is NULL\n")
  cat("rownames(colData) is NULL:", is.null(rownames(colData(l$sce))), "\n")
  cat("colnames(SCE) is NULL:", is.null(colnames(l$sce)), "\n")
}

# Check 2: colnames(SCE) must equal SampleID
cat("\n=== Check 2: colnames(SCE) == SampleID ===\n")
if (!is.null(colnames(l$sce)) && "SampleID" %in% colnames(colData(l$sce))) {
  colnames_sce <- colnames(l$sce)
  sample_ids <- colData(l$sce)$SampleID
  
  cat("colnames(SCE) length:", length(colnames_sce), "\n")
  cat("SampleID length:", length(sample_ids), "\n")
  cat("Identical order:", identical(colnames_sce, sample_ids), "\n")
  cat("Same set:", setequal(colnames_sce, sample_ids), "\n")
  
  if (!identical(colnames_sce, sample_ids)) {
    cat("WARNING: colnames(SCE) != SampleID\n")
    cat("First 5 mismatches:\n")
    mismatches <- data.frame(
      colnames_sce = head(colnames_sce, 5),
      sample_ids = head(sample_ids, 5)
    )
    print(mismatches)
  } else {
    cat("✓ colnames(SCE) == SampleID - PASSED\n")
  }
} else {
  cat("ERROR: Either colnames(SCE) is NULL or SampleID column missing\n")
  cat("colnames(SCE) is NULL:", is.null(colnames(l$sce)), "\n")
  cat("SampleID column exists:", "SampleID" %in% colnames(colData(l$sce)), "\n")
}

# Check 3: SampleID presence and uniqueness
cat("\n=== Check 3: SampleID Quality ===\n")
if ("SampleID" %in% colnames(colData(l$sce))) {
  sample_ids <- colData(l$sce)$SampleID
  
  cat("SampleID column class:", class(sample_ids), "\n")
  cat("Number of unique SampleIDs:", length(unique(sample_ids)), "\n")
  cat("Total number of samples:", length(sample_ids), "\n")
  cat("Any duplicate SampleIDs:", any(duplicated(sample_ids)), "\n")
  cat("Any NA SampleIDs:", any(is.na(sample_ids)), "\n")
  cat("Any empty SampleIDs:", any(sample_ids == ""), "\n")
  
  if (any(duplicated(sample_ids))) {
    cat("WARNING: Duplicate SampleIDs found:\n")
    duplicates <- sample_ids[duplicated(sample_ids)]
    print(table(duplicates))
  }
  
  if (any(is.na(sample_ids))) {
    cat("WARNING: NA SampleIDs found:", sum(is.na(sample_ids)), "\n")
  }
  
  if (any(sample_ids == "")) {
    cat("WARNING: Empty SampleIDs found:", sum(sample_ids == ""), "\n")
  }
} else {
  cat("ERROR: SampleID column not found in colData\n")
}

# Check 4: Test SCE$rn assignment (critical for downstream analysis)
cat("\n=== Check 4: SCE$rn Assignment Test ===\n")
tryCatch({
  # Test if we can assign to SCE$rn
  test_rn <- paste0("test_", seq_len(ncol(l$sce)))
  l$sce$rn <- test_rn
  
  cat("✓ SCE$rn assignment successful\n")
  cat("Assigned rn length:", length(l$sce$rn), "\n")
  cat("First few rn values:", paste(head(l$sce$rn), collapse=", "), "\n")
  
  # Clean up test assignment
  l$sce$rn <- NULL
  cat("✓ Test rn assignment cleaned up\n")
}, error = function(e) {
  cat("ERROR in SCE$rn assignment:", conditionMessage(e), "\n")
  cat("This may cause issues in downstream analysis\n")
})

# Check 5: Comprehensive NULL error prevention
cat("\n=== Check 5: NULL Error Prevention ===\n")
null_checks <- list(
  sce_null = is.null(l$sce),
  colnames_null = is.null(colnames(l$sce)),
  rownames_null = is.null(rownames(l$sce)),
  colData_null = is.null(colData(l$sce)),
  assays_null = is.null(assays(l$sce))
)

cat("NULL checks summary:\n")
for (check_name in names(null_checks)) {
  cat(sprintf("  %s: %s\n", check_name, null_checks[[check_name]]))
}

if (any(unlist(null_checks))) {
  cat("WARNING: NULL values detected that may cause errors\n")
} else {
  cat("✓ No NULL values detected in critical SCE components\n")
}
```

## Expected Output Structure

```
output_directory/
├── annotated_sce.RDS              # Annotated SCE object
├── metadata_sanity_check.xlsx     # Excel summary of annotations
└── validation_log.txt             # Validation and diagnostic output
```

## Customization

### Adding New Experimental Groups
To add new experimental groups, modify the `experimental_groups` list:

```r
experimental_groups <- list(
  # ... existing groups ...
  new_group = list(
    accession = "SH-XX-XXXXXXX X-X",
    timepoint = "Your_timepoint",
    subject = "Your_subject",
    days_since_infusion = your_days,
    reinfusion = TRUE/FALSE
  )
)
```

### Modifying Sentinel Columns
Change the columns used for validation:

```r
# Example: Different sentinel columns
sentinel_cols <- c("roi", "GenesDetected", "Area")

# Example: Single sentinel column
sentinel_cols <- c("roi")
```

### Adding New Metadata Fields
To add new metadata fields, modify the annotation application:

```r
cd[accession == "your_accession",
   `:=`(timepoint = "your_timepoint",
        subject = "your_subject",
        new_field = "new_value")]
```

## Validation and Quality Control

### Required Checks
1. **Row order preservation**: Ensures data integrity during annotation
2. **Sentinel value consistency**: Verifies no unintended changes to key columns
3. **Metadata completeness**: Checks that all samples have required annotations
4. **Object identity**: Verifies raw and annotated objects have same core data
5. **Assay integrity**: Ensures expression data remains unchanged
6. **Column structure**: Validates required metadata columns are present

### Optional Checks
```r
# Check for missing annotations
missing_annotations <- cd[is.na(timepoint) | is.na(subject), .N]
if (missing_annotations > 0) {
  warning(paste("Found", missing_annotations, "samples with missing annotations"))
}

# Check annotation distribution
annotation_summary <- cd[, .N, by = .(timepoint, subject)]
print(annotation_summary)

# Check for duplicate accessions
duplicate_accessions <- cd[duplicated(accession), .N]
if (duplicate_accessions > 0) {
  warning(paste("Found", duplicate_accessions, "duplicate accessions"))
}
```

## Troubleshooting

### Common Issues

1. **Accession not found**: Verify accession identifiers match exactly
2. **Sentinel validation failure**: Check for unintended data modifications
3. **Memory issues**: Process large datasets in chunks if needed
4. **Path errors**: Ensure all file paths are accessible
5. **Column name mismatches**: Check for .dcc suffixes in column names
6. **Assay data corruption**: Verify assay matrices are identical

### Error Handling

The script includes comprehensive validation checks for:
- Data integrity during annotation
- Sentinel value preservation
- Row order consistency
- Object identity verification
- Assay data integrity
- Metadata column structure

## Notes

- Always validate your annotations before proceeding with analysis
- Keep backup copies of raw SCE objects
- Document any manual corrections or special cases
- The Excel output provides a human-readable summary for verification
- The diagnostic output helps identify structural issues with SCE objects
- Validation checks ensure data integrity throughout the annotation process 