# Docker Setup and File Transfer Template - Rech Lab

This template provides a generalized framework for setting up Docker containers and transferring files for GeoMx spatial transcriptomics analysis.

## Overview

This pipeline covers:
- Docker container setup
- File transfer protocols
- Environment configuration
- New instance initialization

## Prerequisites

- Docker installed and configured
- Access to required data files
- Appropriate permissions for file mounting

## Configuration Parameters

### Docker Configuration
```bash
# Container settings
CONTAINER_NAME="your_analysis_container"
IMAGE_NAME="your_docker_image:tag"
WORKING_DIR="/data"
MOUNT_SOURCE="/path/to/your/local/data"
MOUNT_TARGET="/data"
```

### File Transfer Settings
```bash
# Source and destination paths
LOCAL_DATA_DIR="/path/to/your/local/data"
CONTAINER_DATA_DIR="/data"
TRANSFER_FILES=(
    "sce_object.RDS"
    "metadata.csv"
    "analysis_scripts/"
)
```

## Docker Setup

### 1. Pull or Build Docker Image
```bash
# Option 1: Pull existing image
docker pull your_docker_image:tag

# Option 2: Build from Dockerfile
docker build -t your_docker_image:tag .
```

### 2. Create and Run Container
```bash
# Create container with volume mounting
docker run -d \
    --name $CONTAINER_NAME \
    -v $MOUNT_SOURCE:$MOUNT_TARGET \
    -p 8787:8787 \
    -p 8888:8888 \
    $IMAGE_NAME

# Alternative: Interactive mode
docker run -it \
    --name $CONTAINER_NAME \
    -v $MOUNT_SOURCE:$MOUNT_TARGET \
    -p 8787:8787 \
    $IMAGE_NAME
```

### 3. Container Management
```bash
# Start existing container
docker start $CONTAINER_NAME

# Stop container
docker stop $CONTAINER_NAME

# Remove container
docker rm $CONTAINER_NAME

# View container logs
docker logs $CONTAINER_NAME
```

## File Transfer Methods

### Method 1: Volume Mounting (Recommended)
```bash
# Mount local directory to container
docker run -v /local/path:/container/path your_image

# Example for analysis data
docker run -v ~/Desktop/analysis_data:/data your_image
```

### Method 2: Docker Copy
```bash
# Copy files to running container
docker cp local_file.txt $CONTAINER_NAME:/data/

# Copy files from container
docker cp $CONTAINER_NAME:/data/results.txt ./
```

### Method 3: SCP/SFTP (for remote containers)
```bash
# Copy files to remote container
scp -P 2222 local_file.txt user@remote_host:/data/

# Copy files from remote container
scp -P 2222 user@remote_host:/data/results.txt ./
```

## Environment Setup

### 1. Container Environment Variables
```bash
# Set environment variables
docker run -e R_LIBS=/usr/local/lib/R/library \
    -e R_ENVIRON=/usr/local/lib/R/etc/Renviron \
    your_image
```

### 2. R Environment Configuration
```r
# Inside container R session
.libPaths(c("/usr/local/lib/R/library", .libPaths()))
Sys.setenv(R_LIBS_USER = "/usr/local/lib/R/library")
```

### 3. Package Installation
```r
# Install required R packages
install.packages(c("readr", "dplyr", "ggplot2", "pheatmap"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "MAST"))
```

## New Instance Setup

### 1. Initialize New Analysis Instance
```bash
#!/bin/bash
# new-instance-setup.sh

# Configuration
PROJECT_NAME="your_project"
DATA_DIR="/data/$PROJECT_NAME"
SCRIPTS_DIR="/scripts/$PROJECT_NAME"

# Create directory structure
mkdir -p $DATA_DIR/{raw,processed,results}
mkdir -p $SCRIPTS_DIR/{00,01,02}

# Set permissions
chmod 755 $DATA_DIR
chmod 755 $SCRIPTS_DIR

echo "New instance setup complete for $PROJECT_NAME"
```

### 2. File Organization
```
project_structure/
├── data/
│   ├── raw/              # Raw input files
│   ├── processed/        # Intermediate files
│   └── results/          # Final outputs
├── scripts/
│   ├── 00/              # Setup and annotation
│   ├── 01/              # Processing
│   └── 02/              # Analysis
└── docs/                # Documentation
```

### 3. Configuration Files
```bash
# config.sh
export PROJECT_NAME="your_project"
export DATA_DIR="/data/$PROJECT_NAME"
export SCRIPTS_DIR="/scripts/$PROJECT_NAME"
export R_LIBS="/usr/local/lib/R/library"
```

## File Transfer Protocols

### Local Development
```bash
# Mount local development directory
docker run -v $(pwd):/workspace your_image

# Sync changes automatically
docker run -v $(pwd):/workspace:delegated your_image
```

### Remote Development
```bash
# SSH into remote container
ssh -p 2222 user@remote_host

# Use rsync for efficient file transfer
rsync -avz --progress local_dir/ user@remote_host:/data/
```

### Cloud Storage Integration
```bash
# Mount cloud storage (AWS S3 example)
docker run -v ~/.aws:/root/.aws \
    -e AWS_PROFILE=default \
    your_image

# Use AWS CLI inside container
aws s3 cp s3://your-bucket/data/ /data/
```

## Troubleshooting

### Common Issues

1. **Permission denied**: Check file permissions and user mapping
```bash
# Fix permissions
chmod -R 755 /data
chown -R 1000:1000 /data  # Adjust user ID as needed
```

2. **Port conflicts**: Use different ports
```bash
docker run -p 8788:8787 your_image  # Use port 8788 instead of 8787
```

3. **Memory issues**: Increase container memory
```bash
docker run --memory=8g your_image
```

4. **Disk space**: Monitor and clean up
```bash
# Check disk usage
docker system df

# Clean up unused resources
docker system prune
```

### Performance Optimization

```bash
# Use delegated mount for better performance
docker run -v /host/path:/container/path:delegated your_image

# Limit CPU and memory
docker run --cpus=4 --memory=8g your_image

# Use tmpfs for temporary files
docker run --tmpfs /tmp:noexec,nosuid,size=2g your_image
```

## Security Considerations

### File Permissions
```bash
# Restrict file access
chmod 600 sensitive_file.txt
chown root:root sensitive_file.txt

# Use read-only mounts where appropriate
docker run -v /data:/data:ro your_image
```

### Network Security
```bash
# Use internal networks
docker network create analysis_network
docker run --network analysis_network your_image

# Restrict port access
docker run -p 127.0.0.1:8787:8787 your_image
```

## Best Practices

1. **Use specific image tags**: Avoid `latest` tag for reproducibility
2. **Document dependencies**: Keep track of package versions
3. **Backup data**: Regular backups of important files
4. **Monitor resources**: Track CPU, memory, and disk usage
5. **Version control**: Use git for script and configuration management

## Notes

- Always test file transfers with small files first
- Keep container images updated for security
- Document any custom configurations
- Use environment variables for configuration
- Implement proper error handling in scripts 