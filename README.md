# JOINTLY

JOINTLY is an R package for interpretable joint clustering of single-cell and single-nucleus RNA-seq, using kernel-based joint non-negative matrix factorization. JOINTLY effectively captures shared cell states, while robustly accounting for dataset-specific states. It works directly with major single-cell genomics frameworks such as [Seurat](https://github.com/satijalab/seurat) or [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment) (or Matrix-like) objects.


# Installation

valiDrops can be installed directly from GitHub using either [devtools](https://cran.r-project.org/web/packages/devtools/index.html) or [remotes](https://cran.r-project.org/web/packages/remotes/index.html). 

```{R}
install.packages("remotes")
remotes::install_github("madsen-lab/rJOINTLY")
```


# Usage

JOINTLY uses a list of scRNA-seq expression objects such as [Seurat](https://github.com/satijalab/seurat), [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment) or raw expression matrices. 
Here we demonstrate how to make a JOINTLY object and how to run the clustering algorithm. 


```{R}
## Load libraries
library(Seurat)
library(rJOINTLY)

## Load test data
seu <- readRDS(file = "../data/data_seurat.rds")

## Run JOINTLY with default parameters
batch_var = 'donor_id'
out <- jointly(seu, batch.var = batch_var)

## Run UMAP
jointly_out <- out[[1]]
jointly_out <- RunUMAP(jointly_out, reduction = "JOINTLY", dims = 1:10)
DimPlot(jointly_out, label = TRUE, group.by = batch_var)
```

### Citation
_Møller AF, Madsen JGS et al. **Interpretable joint clustering of single-cell transcriptomes (2023)** Unpublished_  <br/>
