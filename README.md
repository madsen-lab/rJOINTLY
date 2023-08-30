# JOINTLY

JOINTLY is an R package for interpretable joint clustering of single-cell and single-nucleus RNA-seq, using kernel-based joint non-negative matrix factorization. JOINTLY effectively captures shared cell states, while robustly accounting for dataset-specific states. It works directly with major single-cell genomics frameworks such as [Seurat](https://github.com/satijalab/seurat) or [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment) (or Matrix-like) objects.

JOINTLY is also available as a Python package [here](https://github.com/madsen-lab/pyJOINTLY)

Scripts for reproducing the analyses in the manuscript are available [here](https://github.com/madsen-lab/JOINTLY_reproducibility)

For the white adipose tissue atlas (WATLAS), the model weights are availiable [here](https://zenodo.org/deposit/8086433) and the atlas can be explored and downloaded [here](https://singlecell.broadinstitute.org/single_cell/study/SCP2289/an-integrated-single-cell-and-single-nucleus-rna-seq-white-adipose-tissue-atlas-watlas)

# Installation

JOINTLY can be installed directly from GitHub using either [devtools](https://cran.r-project.org/web/packages/devtools/index.html) or [remotes](https://cran.r-project.org/web/packages/remotes/index.html). 

```{R}
install.packages("remotes")
remotes::install_github("madsen-lab/rJOINTLY")
```

JOINTLY depends on R version >= 4.1.0 has been tested on Unix, Windows and macOS systems. It has several package dependencies, which should be installed automatically. In case they are not, the dependencies are listed here: BiocParallel, Matrix, R.utils, Rcpp, Seurat, SharedObject, SingleCellExperiment, e1071, inflection, irlba, scran, scry, and transformGamPoi.

# Basic usage

JOINTLY uses a list of scRNA-seq expression objects such as [Seurat](https://github.com/satijalab/seurat), [SingleCellExperiment](https://github.com/drisso/SingleCellExperiment) or raw expression matrices. 
Here, we demonstrate how to obtain a testing dataset and how to perform joint clustering and interpretation using JOINTLY. In the example below, we use SeuratData to download a test dataset, Seurat to handle the dataset, aricode to calculate ARI between cell type labels and cluster labels, and UCell to calculate module scores:

```{R}
### How to setup the environment for testing JOINTLY
## Install libraries from CRAN
install.packages("remotes")
install.packages("Seurat")
install.packages("aricode")

## Install libraries from GitHub
remotes::install_github("madsen-lab/rJOINTLY")
remotes::install_github("carmonalab/UCell")

## Load libraries
library(Seurat)
library(aricode)
library(UCell)
library(JOINTLY)

## Load the liver dataset
# Datasets, labels, embeddings and results used for testing JOINTLY can be downloaded here: https://zenodo.org/record/8298157
liver <- readRDS("Human_Liver.rds")

### Joint clustering and evaluation with JOINTLY
## Run JOINTLY
# To run in parallel pass the option: bpparam = BiocParallel::MulticoreParam()
# The resulting object is a list that holds a Seurat object as the first entry and module as the second entry
liver <- jointly(liver, batch.var = "batch_label")

## Import labels and add to the object
labels <- read.delim("Human_Liver_labels.tsv")
labels <- labels[ match(colnames(liver$object), labels$X),]
liver$object$celltype <- labels$Transferred_labels

## Subset the object to only named cell types with more than 10 cells
liver$object <- subset(liver$object, celltype %in% names(which(table(liver$object$celltype) >= 10)))
liver$object <- subset(liver$object, celltype != "")

## Run clustering
liver$object <- FindNeighbors(liver$object, reduction = "JOINTLY", dims = 1:15)
liver$object <- FindClusters(liver$object, res = 0.05)

## Calculate ARI and NMI
aricode::ARI(liver$object$celltype, liver$object$seurat_clusters) # 0.7415981
aricode::NMI(liver$object$celltype, liver$object$seurat_clusters) # 0.6325386

## Embed using UMAP
liver$object <- RunUMAP(liver$object, reduction = "JOINTLY", dims = 1:15)

## Plot by cluster, cell type and batch
# Setting a seed and re-running produces a different result. Sometimes it is worth running JOINTLY up to 5 times and choosing the best solution.
DimPlot(liver$object, group.by = "seurat_clusters", label = TRUE) + DimPlot(liver$object, group.by = "celltype", label = TRUE) + DimPlot(liver$object, group.by = "batch_label", label = FALSE) & NoLegend()

### Interpretation of JOINTLY
## Calculate module scores
liver$object <- AddModuleScore_UCell(liver$object, features = liver$modules)

## Visualize module scores on the UMAP
FeaturePlot(liver$object, features = "factor_1_UCell")
```

# Advanced usage

It is possible to use the learned H and F matrices in JOINTLY to perform batch-aware smoothening of expression values using an approach similar to [MAGIC](https://pubmed.ncbi.nlm.nih.gov/29961576/). Please be aware, that we advise to use these values with caution, and only for visualization purposes. In order to do so, it is necessary to run JOINTLY without using the wrapper function:

```{R}
## Install and load required libraries
remotes::install_github("immunogenomics/presto")
install.packages("expm")
library("expm")
library("presto")

## Reload the liver dataset
liver <- readRDS("Human_Liver.rds")
liver <- readRDS("C:/Users/jgsm/OneDrive - Syddansk Universitet/Human_Liver.rds")

## Run JOINTLY without using the wrapper function
proc <- preprocess(liver, "batch_label")
cpca <- cpca(proc)
inputs <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca)

## Insert the H matrix as a DimReduc object in the original Seurat object
Hmat <- solved$Hmat.scaled
Hmat <- Hmat[ match(colnames(liver), rownames(Hmat)),]
liver[["JOINTLY"]] <- Seurat::CreateDimReducObject(Hmat, assay = "RNA")

## Embed using UMAP
liver <- RunUMAP(liver, reduction = "JOINTLY", dims = 1:15)

## Calculate smoothend expression values for each dataset
## To allow more T steps change the n_steps variable
n_steps = 1
recon_expression_list <- list()
for (idx in 1:length(proc)){
  recon_kernel = solved$Fmat[[idx]] %*% solved$Hmat[[idx]]
  recon_kernel = 0.5 * (recon_kernel + t(recon_kernel))
  D = 1 / rowSums(recon_kernel)
  T = matrix(diag(D),ncol=length(D)) %*% recon_kernel
  T_steps = T %^% n_steps
  recon_expression = as.data.frame(T_steps %*% t(cpca$normalized[[idx]]))
  rownames(recon_expression) = colnames(cpca$normalized[[idx]])
  recon_expression_list[[idx]] = recon_expression
}

## Combine across batches
recon_expression = as.data.frame(do.call("rbind", recon_expression_list))

## Create a new assay with the imputed expression values
liver[["imputed"]] <- CreateAssayObject(data = t(recon_expression[match(colnames(liver), rownames(recon_expression)),]))

## Visualize IGF1 expression on the UMAP using raw and imputed values
# IGF1 is primarily expressed by hepatocytes
fp_raw <- FeaturePlot(liver, "ENSG00000017427", reduction = "umap")
DefaultAssay(liver) <- "imputed"
fp_imputed <- FeaturePlot(liver, "ENSG00000017427", reduction = "umap")
fp_raw + fp_imputed

## Import labels and add to the object
labels <- read.delim("Human_Liver_labels.tsv")
labels <- labels[ match(colnames(liver), labels$X),]
liver$celltype <- labels$Transferred_labels

## Calculate pseudo-bulk expression in hepatocytes for each donor
liver$group <- paste(liver$celltype, liver$batch_label, sep="_")
psb_imputed <- AverageExpression(liver, assays = "imputed", group.by = "group")$imputed
psb_raw <- AverageExpression(liver, assays = "RNA", group.by = "group")$RNA

## Calculate coefficient of variation for IGF1 as an example
# Batch-aware imputed values have a much smaller coefficient of variation indicating batch effects have been reduced
igf1_hep_imputed <- psb_imputed[ rownames(psb_imputed) == "ENSG00000017427", grep("Hepatocytes", colnames(psb_imputed))]
igf1_hep_raw <- psb_raw[ rownames(psb_raw) == "ENSG00000017427", grep("Hepatocytes", colnames(psb_raw))]
sd(igf1_hep_raw) / mean(igf1_hep_raw) # 1.045458
sd(igf1_hep_imputed) / mean(igf1_hep_imputed) # 0.5184951

## Import labels and subset the object
labels <- read.delim("Human_Liver_labels.tsv")
labels <- labels[ match(colnames(liver$object), labels$X),]
liver$object$celltype <- labels$Transferred_labels
liver$object <- subset(liver$object, celltype %in% names(which(table(liver$object$celltype) >= 10)))
liver$object <- subset(liver$object, celltype != "")

## Calculate coefficient of variation for all genes in all cell types
cov.raw <- c()
cov.imputed <- c()
for (ct in unique(liver$celltype)) {
  imputed <- psb_imputed[ , grep(ct, colnames(psb_imputed))]
  cov.imputed <- c(cov.imputed, apply(imputed,1,FUN="sd") / apply(imputed,1,FUN="mean"))
  raw <- psb_raw[, grep(ct, colnames(psb_raw))]
  cov.raw <- c(cov.raw, apply(raw,1,FUN="sd") / apply(raw,1,FUN="mean"))
}
median(cov.raw, na.rm=TRUE) # 1.38173
median(cov.imputed, na.rm=TRUE) # 0.513786

## Find marker genes using presto using raw nad imputed values
res.raw <- wilcoxauc(liver, group_by = "celltype", seurat_assay = "RNA")
res.raw$id <- paste(res.raw$group, res.raw$feature, sep="_")
res.imputed <- wilcoxauc(liver, group_by = "celltype", seurat_assay = "imputed")
res.imputed$id <- paste(res.imputed$group, res.imputed$feature, sep="_")

## Define list of markers per cell type found using the raw values
markers.raw <- res.raw[ res.raw$logFC >= log2(1.1) & res.raw$padj <= 0.05 & res.raw$auc >= 0.6,"id"]

## Compare the average AUC for the marker genes between the raw and imputed values
# The average AUC for marker genes found in the raw data is much higher after imputation indicating better cell type seperation.
mean(res.raw[ res.raw$id %in% markers.raw,"auc"]) # 0.6866768
mean(res.imputed[ res.imputed$id %in% markers.raw,"auc"]) # 0.8296886
```

# Citation
If you use JOINTLY in your work, please consider citing our manuscript:

_Møller AF, Madsen JGS et al. **Interpretable joint clustering of single-cell transcriptomes (2023)** Unpublished_  <br/>
