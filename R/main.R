#' @title jointly
#'
#' @description Function to cluster single-cell and single-nucleus RNA-seq data using JOINTLY. 
#' 
#' @param data A single (or list of) Seurat, SingleCellExperiment, or Matrix-like object(s) containing raw counts.
#' @param batch.var A variable to split Seurat or SingleCellExperiment objects into a list [default = NULL]
#' @param factors The number of factors to calculate [default = 20]
#' @param nfeat The number of features to select [default = 1000]
#' @param modules Boolean (TRUE or FALSE) determining whether or not to return gene modules per factor [default = TRUE].
#' @param init The method for initializing the H and F matrices ("random" or "clustering").[default = "clustering"]
#' @param selection.method The method to use selecting highly-variable features [default = "deviance"]
#' @param decay.k The number of neighbors to use for the decay function [default = 5]
#' @param decay.alpha The power of the decay function [default = 5] 
#' @param cpca.threshold The minimum amount of variance to be explained by common PCs before adding individual PCs [default = 0.8]
#' @param cpca.kc The number of common PCs to calculate [default = 20]
#' @param cpca.ki The number of invididual PCs to calcuate [default = 20]
#' @param alpha.loss Alpha parameter of the loss function [default = 1]
#' @param mu.loss Mu parameter of the loss function [default = 20]
#' @param lambda.loss Lambda parameter of the loss function [default = 1]
#' @param beta.loss Beta parameter of the loss function [default = 5]
#' @param snn.k The number of neighbors to use for building the SNN graph [default = 30]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#' @param iter.max The maximum number of iterations to perform [default = 100]
#' @param verbose Boolean (TRUE or FALSE) determining verbosity [default = TRUE]
#' @param bpparam *Param to use for parallel processing [default = SerialParam()]
#' @param ... Additional parameters to pass to functions within JOINTLY
#'
#' @return Returns either a Seurat object, a SingleCellExperiment object or a list containing the H and W matrices
#' @export
#' @import R.utils
#' @import Rcpp
#' @import inflection
#' @importFrom methods as
#' @importFrom stats coef lm.fit rnorm runif
#' @useDynLib JOINTLY
#'

jointly <- function(data, batch.var = NULL, factors = 20, nfeat = 1000, modules = TRUE, init = "clustering", bpparam = SerialParam(), selection.method = "deviance", decay.k = 5, decay.alpha = 5, cpca.threshold = 0.8, cpca.kc = 20, cpca.ki = 20, alpha.loss = 1, mu.loss = 20, lambda.loss = 1, beta.loss = 5, snn.k = 30, ncpu = 1, iter.max = 100, verbose = TRUE, ...) {
  # TODO: Check parameters
  # TODO: Check for duplicated barcode names
  # TODO: Check for missing names for datatsets already provided as a list
  
  # Preprocess
  preprocessed <- R.utils::doCall(JOINTLY::preprocess, args = ..., alwaysArgs = list(data = data, batch.var = batch.var, verbose = verbose))
  
  # CPCA
  cpca.res <- R.utils::doCall(JOINTLY::cpca, args = ..., alwaysArgs = list(bpparam = bpparam, dataset.list = preprocessed, nfeat = nfeat, selection.method = selection.method, threshold = cpca.threshold, kc = cpca.kc, ki = cpca.ki, ncpu = ncpu, iter.max = iter.max, verbose = verbose))
  norm.list <- cpca.res$normalized
  cpca.list <- cpca.res$cpca 
  
  # prepareData
  inputs <- R.utils::doCall(JOINTLY::prepareData, args = ..., alwaysArgs = list(dataset.list = cpca.list, k.decay = decay.k, alpha = decay.alpha, k.snn = snn.k, verbose = verbose))
  kernel.list <- inputs$kernels
  snn.list <- inputs$snn
  rare.list <- inputs$rareity
  
  # Solve
  mat <- R.utils::doCall(JOINTLY::JOINTLYsolve, args = ..., alwaysArgs = list(bpparam = bpparam, kernel.list = kernel.list, snn.list = snn.list, rare.list = rare.list, cpca.result = cpca.res, k = factors, init = init, iter.max = iter.max, alpha = alpha.loss, mu = mu.loss, lambda = lambda.loss, beta = beta.loss, ncpu = ncpu, progressbar = verbose))
  
  ## Finalize results
  Hmat <- mat$Hmat.scaled
  if (is.null(batch.var)) {
    object.list <- list()
    for (ds in 1:length(data)) {
      tmp <- data[[ds]]
      Hmat.tmp <- Hmat[ rownames(Hmat) %in% colnames(data[[ds]]),]
      Hmat.tmp <- Hmat.tmp[ match(colnames(data[[ds]]), rownames(Hmat.tmp)),]
      if (class(tmp)[1] == "Seurat") {
        tmp[["JOINTLY"]] <- Seurat::CreateDimReducObject(Hmat.tmp, assay = "RNA")
      } else if (class(tmp)[1] == "SingleCellExperiment") {
        SingleCellExperiment::reducedDim(tmp, "JOINTLY") <- SingleCellExperiment::reduced.dim.matrix(Hmat.tmp)
      } else {
        tmp <- Hmat.tmp
      }
      object.list[[ds]] <- tmp
    }
  } else {
    Hmat <- Hmat[ match(colnames(data), rownames(Hmat)),]
    if (class(data)[1] == "Seurat") {
      data[["JOINTLY"]] <- Seurat::CreateDimReducObject(Hmat, assay = "RNA")
      object.list <- data
    } else if (class(data)[1] == "SingleCellExperiment") {
      SingleCellExperiment::reducedDim(data, "JOINTLY") <- SingleCellExperiment::reduced.dim.matrix(Hmat)
      object.list <- data
    } else {
      object.list <- Hmat
    }
  }
  
  ## Process modules
  if (modules) {
    Wmat <- mat$Wmat.scaled
    modules <- list()
    for (i in 1:factors) {
      modules[[length(modules)+1]] <- names(sort(Wmat[,i], decreasing = TRUE))[1:inflection::uik(y = sort(Wmat[,i], decreasing = TRUE), x = seq(1, nrow(Wmat),1))]
      names(modules)[length(modules)] <- paste("factor_", i, sep="")
    }
    result.list = list(object = object.list, modules = modules)
  } else {
    result.list = list(object = object.list)
  }
  
  # Return
  return(result.list)
}
