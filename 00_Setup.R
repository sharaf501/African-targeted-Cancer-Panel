# 00_setup.R
# Install and load required packages

# CRAN packages ----
cran_packages <- c(
  "tidyverse",      # Data manipulation
  "data.table",     # Fast data reading
  "rentrez",        # PubMed API access
  "europepmc",      # Europe PMC API
  "httr",           # HTTP requests
  "jsonlite",       # JSON parsing
  "xml2",           # XML parsing
  "rvest",          # Web scraping
  "BiocManager",    # Bioconductor installer
  "ggplot2",        # Visualization
  "plotly",         # Interactive plots
  "pheatmap",       # Heatmaps
  "ComplexHeatmap", # Advanced heatmaps
  "VennDiagram",    # Venn diagrams
  "UpSetR",         # UpSet plots
  "gridExtra",      # Multiple plots
  "knitr",          # Reports
  "rmarkdown",      # Markdown reports
  "DT",             # Interactive tables
  "writexl",        # Excel export
  "openxlsx",       # Excel manipulation
  "scales",         # Scale functions
  "RColorBrewer",   # Color palettes
  "viridis",        # Color palettes
  "cowplot",        # Plot arrangements
  "patchwork",      # Combine plots
  "ggsci",          # Scientific color palettes
  "ggrepel",        # Label repelling
  "corrplot",       # Correlation plots
  "reshape2",       # Data reshaping
  "stringr",        # String manipulation
  "lubridate",      # Date handling
  "forcats",        # Factor handling
  "glue",           # String interpolation
  "here",           # Path management
  "fs"              # File system operations
)

# Install CRAN packages ----
install.packages(cran_packages[!cran_packages %in% installed.packages()])

# Bioconductor packages
bioc_packages <- c(
  "Biostrings",           # Sequence analysis
  "GenomicRanges",        # Genomic intervals
  "VariantAnnotation",    # Variant handling
  "biomaRt",              # Ensembl access
  "maftools",             # MAF file analysis
  "TCGAbiolinks",         # TCGA data
  "org.Hs.eg.db",         # Gene annotation
  "GO.db",                # Gene Ontology
  "KEGG.db",              # KEGG pathways
  "reactome.db",          # Reactome pathways
  "AnnotationDbi",        # Annotation queries
  "BSgenome.Hsapiens.UCSC.hg38", # Reference genome
  "TxDb.Hsapiens.UCSC.hg38.knownGene", # Gene models
  "ensembldb",            # Ensembl databases
  "oncoPredict",          # Drug response prediction
  "clusterProfiler",      # Enrichment analysis
  "DOSE",                 # Disease ontology
  "enrichplot"            # Enrichment visualization
)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioc_packages[!bioc_packages %in% installed.packages()])

# Load core packages ----
library(tidyverse)
library(data.table)
library(here)
library(fs)
library(tinytex)

# Create project structure ----
dir_create(here("data", "raw"))
dir_create(here("data", "processed"))
dir_create(here("data", "external"))
dir_create(here("scripts"))
dir_create(here("results", "figures"))
dir_create(here("results", "tables"))
dir_create(here("results", "reports"))
dir_create(here("functions"))
dir_create(here("reports"))

# Set global options ----
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  timeout = 300
)

# Create configuration file ----
config <- list(
  project_name = "African targeted PanCancer Gene Panel",
  author = "Sharafudeen Abubakar",
  date = Sys.Date(),
  ncbi_api_key = "", # Get from https://www.ncbi.nlm.nih.gov/account/settings/
  use_cbioportal = TRUE,
  use_civic = TRUE,
  use_icgc = TRUE,
  use_gnomad = TRUE
)

saveRDS(config, here("data", "config.rds"))
