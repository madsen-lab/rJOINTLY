# JOINTLY

JOINTLY is an R package for interpretable joint clustering of single-cell and single-nucleus RNA-seq, using kernel-based joint non-negative matrix factorization. JOINTLY effectively captures shared cell states, while robustly accounting for dataset-specific states. It works directly on Seurat or SingleCellExperiment (or Matrix-like) objects.


# Installation

valiDrops can be installed directly from GitHub using either `[devtools](https://cran.r-project.org/web/packages/devtools/index.html)` or `[remotes](https://cran.r-project.org/web/packages/remotes/index.html)`. 

```{R}
install.packages("remotes")
remotes::install_github("madsen-lab/rJOINTLY")
```


# Usage

JOINTLY uses a list of scRNA-seq expression objects e.g. [scanpy](https://github.com/scverse/scanpy) or raw expression matrices as numpy arrays. 
Here we demonstrate how to make a JOINTLY object and how to run the clustering algorithm. 


```{R}
## Load libraries
import scanpy as sc
from jointly import jointly, CreateJointlyObject_from_scanpyList


## Load test data
data_list = [sc.read_h5ad('pancdata/panc1.h5ad'), 
             sc.read_h5ad('pancdata/panc2.h5ad'), 
             sc.read_h5ad('pancdata/panc3.h5ad')]


## Preprocess Seurat a object
proc <- preprocess(seu, batch_var)

## Run cPCA
cpca <- cpca(proc)


## Prepare kernel and SNN graph
inputs <- prepareData(cpca$cpca)

## Running JOINTLY with default parameters
jointly(jointlyObject)

## Running JOINTLY with default parameters
solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), progressbar = FALSE)

## Get and scale clustering matrix
H <- t(do.call("cbind",solved$Hmat))
H <- H[ match(colnames(data), rownames(H)),]
H <- scale(H)
H <- t(scale(t(H)))
colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")

## Cluster the matrix using e.g. hierarchical clustering
N_clusters = 5
snn <- SNN.Construction(H, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = N_clusters)
```

> _For any requests, please reach out to: <br/>[Jesper Madsen](jgsm@imada.sdu.dk) or [Andreas Møller](andreasfm@bmb.sdu.dk)_


### Citation
_Møller AF, Madsen JGS et al **Interpretable joint clustering of single-cell transcriptomes (2023)** Unpublished_  <br/>
