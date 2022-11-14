#' @title jointly
#'
#' @description Function to cluster single-cell and single-nucleus RNA-seq data using JOINTLY. 
#' 
#' @param data A single (or list of) Seurat, SingleCellExperiment, or Matrix-like object(s) containing raw counts.
#' @param batch.var A variable to split Seurat or SingleCellExperiment objects into a list [default = NULL]
#' @param factors The number of factors to calculate [default = 20]
#' @param nfeat The number of features to select [default = 1000]
#' @param init The method to use for initializing H and F matrices ("random" or "clustering"). [default = "clustering"]
#' @param selection.method The method to use selecting highly-variable features [default = "deviance"]
#' @param decay.k The number of neighbors to use for the decay function [default = 5]
#' @param decay.alpha The power of the decay function [default = 1] 
#' @param cpca.threshold The minimum amount of variance to be explained by common PCs before adding individual PCs [default = 0.8]
#' @param cpca.kc The number of common PCs to calculate [default = 20]
#' @param cpca.ki The number of invididual PCs to calcuate [default = 20]
#' @param alpha.loss Alpha parameter of the loss function [default = 100]
#' @param mu.loss Mu parameter of the loss function [default = 1]
#' @param lambda.loss Lambda parameter of the loss function [default = 100]
#' @param beta.loss Beta parameter of the loss function [default = 1]
#' @param snn.k The number of neighbors to use for building the SNN graph [default = 50]
#' @param rare.k The number of neighbors to calculate rareity scores [default = 10]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#' @param iter.max The maximum number of iterations to perform [default = 100]
#' @param verbose Boolean (TRUE or FALSE) determining verbosity [default = TRUE]
#' @param ... Additional parameters to pass to functions within JOINTLY
#'
#' @return The same type of object as inputted to the function with added an added dimensional reduction called JOINTLY, which can be used for embedding and clustering.
#' @export
#' @import R.utils
#' @import Rcpp
#' @importFrom methods as
#' @importFrom stats coef lm.fit rnorm runif
#' @useDynLib JOINTLY
#'

jointly <- function(data, batch.var = NULL, factors = 20, nfeat = 1000, init = "clustering", selection.method = "deviance", decay.k = 5, decay.alpha = 2, cpca.threshold = 0.8, cpca.kc = 20, cpca.ki = 20, alpha.loss = 10, mu.loss = 10, lambda.loss = 5, beta.loss = 10, snn.k = 100, rare.k = 5, ncpu = 1, iter.max = 100, verbose = TRUE, ...) {
  # TODO: Check parameters
  
  ## Preprocess
  if (verbose) { message("Preprocessing dataset.")}
  preprocessed <- R.utils::doCall(JOINTLY::preprocess, args = ..., alwaysArgs = list(data = data, batch.var = batch.var))
  
  ## CPCA
  if (verbose) { message("Computing consensus PCA.")}
  cpca.res <- R.utils::doCall(JOINTLY::cpca, args = ..., alwaysArgs = list(dataset.list = preprocessed, nfeat = nfeat, selection.method = selection.method, threshold = cpca.threshold, kc = cpca.kc, ki = cpca.ki, ncpu = ncpu, iter.max = iter.max, verbose = verbose))
  norm.list <- cpca.res$normalized
  cpca.list <- cpca.res$cpca 
  
  ## prepareData
  if (verbose) { message("Computing decay kernels, SNN graphs and rareity scores.")}
  inputs <- R.utils::doCall(JOINTLY::prepareData, args = ..., alwaysArgs = list(dataset.list = cpca.list, k.decay = decay.k, alpha = decay.alpha, k.rare = rare.k, k.snn = snn.k))
  kernel.list <- inputs$kernels
  snn.list <- inputs$snn
  rare.list <- inputs$rareity
  
  ## Solve
  if (verbose) { message("Solving matrices.")}
  mat <- R.utils::doCall(JOINTLY::JOINTLYsolve, args = ..., alwaysArgs = list(kernel.list = kernel.list, snn.list = snn.list, rare.list = rare.list, norm.list = norm.list, k = factors, iter.max = iter.max, alpha = alpha.loss, mu = mu.loss, lambda = lambda.loss, beta = beta.loss, ncpu = ncpu))
  
  ## Finalize results
  # Scale the H matrices
  if (verbose) { message("Finalizing.")}
  Hmat <- mat$Hmat
  for (i in 1:length(Hmat)) {	
    Hmat[[i]] <- t(scale(t(Hmat[[i]])))
    rownames(Hmat[[i]]) <- paste("JOINTLY_", seq(1,factors,1), sep="")
  }
  
  # Insert results into input objects
  if (is.null(batch.var)) {
    result.list <- list()
    names(result.list) <- names(data)
    for (ds in 1:length(data)) {
      tmp <- data[[ds]]
      res <- t(Hmat[[which(names(Hmat) == names(data)[[ds]])]])
      res <- res[ match(colnames(tmp), rownames(res)),]
      colnames(res) <- paste("JOINTLY_", seq(1, ncol(res),1), sep="")
      if (class(tmp)[1] == "Seurat") {
        tmp[["JOINTLY"]] <- Seurat::CreateDimReducObject(res, assay = "RNA")
      } else if (class(tmp)[1] == "SingleCellExperiment") {
        SingleCellExperiment::reducedDim(tmp, "JOINTLY") <- SingleCellExperiment::reduced.dim.matrix(res)
      } else {
        tmp <- res
      }
      result.list[[ds]] <- tmp
    }
  } else {
    res <- t(do.call("cbind", Hmat))
    res <- res[ match(colnames(data), rownames(res)),]
    colnames(res) <- paste("JOINTLY_", seq(1, ncol(res),1), sep="")
    if (class(data)[1] == "Seurat") {
      data[["JOINTLY"]] <- Seurat::CreateDimReducObject(res, assay = "RNA")
      result.list <- data
    } else if (class(data)[1] == "SingleCellExperiment") {
      SingleCellExperiment::reducedDim(data, "JOINTLY") <- SingleCellExperiment::reduced.dim.matrix(res)
      result.list <- data
    } else {
      result.list <- res
    }
  }
  
  # Return
  return(result.list)
}
