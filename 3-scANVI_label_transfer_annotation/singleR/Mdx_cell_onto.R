# Install and load BiocManager if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# Install necessary Bioconductor packages
#BiocManager::install(c("SingleR", "scRNAseq", "SummarizedExperiment", "MatrixGenerics", "GenomicRanges", "GenomeInfoDb"))
#py_install(c("rpy2", "anndata2ri"), pip = TRUE)

# Load necessary libraries
library(SingleR)
library(scRNAseq)
library(SummarizedExperiment)
library(MatrixGenerics)
library(reticulate)
library(zellkonverter)
library(ontoProc)


Sys.unsetenv("RETICULATE_PYTHON")

#reticulate::use_python("/Users/srivalli/conda/envs/scanpy_new/bin/python")
#py_config()
use_condaenv("/Users/srivalli/conda/envs/scanpy_new", required = TRUE)
py_config()

anndata2ri <- import("anndata2ri")
rpy2 <- import("rpy2")

anndata2ri$activate()
pandas2ri <- import("rpy2.robjects.pandas2ri")

# Load your h5ad file using Scanpy in Python
sc <- import('scanpy')
adata <- readH5AD('/Users/srivalli/Documents/GitHub/hofmann_dmd/data/heart_mm_nuclei_CB_QC_23_05_2024.raw.h5ad')
mdx_data <- adata[,adata$genotype == "Mdx"]
dim(mdx_data)

# Log-transform the normalized counts
normalized_counts <- assay(mdx_data, "sqrt_norm") 

# Perform log transformation
logcounts <- log1p(normalized_counts)

# Add the log-transformed counts as a new assay
assay(mdx_data, "logcounts") <- logcounts

# Assuming 'adata' is your SingleCellExperiment object
ref <-  celldex::MouseRNAseqData()

#Traceback
cl <- getOnto("cellOnto")
cl

#Translate ontologies
translated <- cl$name[ref$label.ont]
head(translated)

# Subset to only include cardiac data
#ref <- ref[, ref$label.fine == "Cardiac"]

pred <- SingleR(test = mdx_data, ref = ref, labels = translated)

# Convert prediction to DataFrame and assign to adata
pred_df <- as.data.frame(pred)
mdx_data$SingleR <- pred_df$labels

# Save annotated data back to h5ad format
writeH5AD(mdx_data,'/Users/srivalli/Documents/GitHub/hofmann_dmd/data/annotated_Mdx_ont_data.h5ad')