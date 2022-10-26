#' @title jointly
#'
#' @description Function to cluster single-cell and single-nucleus RNA-seq data using JOINTLY. 
#' 
#' @param data A single (or list of) Seurat, SingleCellExperiment, or Matrix-like object(s) containing raw counts.
#' @param batch.var A variable to split Seurat or SingleCellExperiment objects into a list [default = NULL]
#' @param factors The number of factors to calculate [default = 20]
#' @param nfeat The number of features to select [default = 1000]
#' @param selection.method The method to use selecting highly-variable features [default = "deviance"]
#' @param decay.k The number of neighbors to use for the decay function [default = 5]
#' @param decay.alpha The power of the decay function [default = 2] 
#' @param cpca.threshold The minimum amount of variance to be explained by common PCs before adding individual PCs [default = 0.8]
#' @param cpca.kc The number of common PCs to calculate [default = 20]
#' @param cpca.ki The number of invididual PCs to calcuate [default = 20]
#' @param alpha.loss Alpha parameter of the loss function [default = 100]
#' @param mu.loss Mu parameter of the loss function [default = 1]
#' @param lambda.loss Lambda parameter of the loss function [default = 100]
#' @param beta.loss Beta parameter of the loss function [default = 1]
#' @param snn.k The number of neighbors to use for building the SNN graph [default = "rice" (= (2 * n^(1/3)))]
#' @param snn.add The number of neighbors to add in addition to k to use for building the SNN graph [default = 20]
#' @param snn.rare A boolean (TRUE or FALSE) indicating if rare cell index should be calculated from the distance to the Kth neighbour during SNN graph construction [default = TRUE]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#' @param iter.max The maximum number of iterations to perform [default = 100]
#' @param verbose Boolean (TRUE or FALSE) determining verbosity [default = TRUE]
#' @param ... Additional parameters to pass to functions within JOINTLY
#'
#' @return A list containing a list (per-batch) of H matrices, a list (per-batch) of F matrices and a list (per-batch) of W matrices.
#' @export
#' @import R.utils
#' @import Rcpp
#' @importFrom methods as
#' @importFrom stats coef lm.fit rnorm runif
#' @useDynLib JOINTLY
#'

jointly <- function(data, batch.var = NULL, factors = 20, nfeat = 1000, selection.method = "deviance", decay.k = 5, decay.alpha = 2, cpca.threshold = 0.8, cpca.kc = 20, cpca.ki = 20, alpha.loss = 100, mu.loss = 1, lambda.loss = 100, beta.loss = 1, snn.k = "rice", snn.add = 20, snn.rare = TRUE, ncpu = 1, iter.max = 100, verbose = TRUE, ...) {
  # TODO: Check parameters
  
  # Preprocess
  if (verbose) { message("Preprocessing dataset.")}
  preprocessed <- R.utils::doCall(JOINTLY::preprocess, args = ..., alwaysArgs = list(data = data, batch.var = batch.var))
  
  # CPCA
  if (verbose) { message("Computing consensus PCA.")}
  cpca.res <- R.utils::doCall(JOINTLY::cpca, args = ..., alwaysArgs = list(dataset.list = preprocessed, nfeat = nfeat, selection.method = selection.method, threshold = cpca.threshold, kc = cpca.kc, ki = cpca.ki, ncpu = ncpu, iter.max = iter.max, verbose = verbose))
  norm.list <- cpca.res$normalized
  cpca.list <- cpca.res$cpca 
  
  # alphaDecay kernel
  if (verbose) { message("Computing alpha decay kernel.")}
  kernel.list <- R.utils::doCall(JOINTLY::alphaDecay, args = ..., alwaysArgs = list(dataset.list = cpca.list, k = decay.k, alpha = decay.alpha))
  
  # SNN graphs
  if (verbose) { message("Computing SNN graph.")}
  snn.res <- R.utils::doCall(JOINTLY::computeSNN, args = ..., alwaysArgs = list(dataset.list = cpca.list, k = snn.k, add = snn.add, rare = snn.rare))
  snn.list <- snn.res$snn
  rare.list <- snn.res$rare
  
  # Solve
  if (verbose) { message("Solving matrices.")}
  mat <- R.utils::doCall(JOINTLY::JOINTLYsolve, args = ..., alwaysArgs = list(kernel.list = kernel.list, snn.list = snn.list, rare.list = rare.list, norm.list = norm.list, k = factors, iter.max = iter.max, alpha = alpha.loss, mu = mu.loss, lambda = lambda.loss, beta = beta.loss, ncpu = ncpu))

  # Return
  return(mat)
}
