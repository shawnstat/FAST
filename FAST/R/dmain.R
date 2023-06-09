#' Flexible Analysis of Spatial Transcriptomics (FAST) Data: A Deconvolution Approach
#'
#' @docType package
#' @name FAST


#' Default configuration for FAST
#'
#' A list with parameters
#'
#' r: integer; rank of decomposition
#'
#' lambda_1: numeric; tuning parameter for the graph matrix
#'
#' lambda_2: numeric; tuning parameter for sparsity levels of W
#'
#' n_iter: integer; number of iterations
#'
#' converge: numeric; convergence criteria
#'
#' @examples
#' # display all default settings
#' dmain_config
#'
#' # create a new settings object with r set to 5
#' new_config <- dmain_config
#' new_config$r <- 5
#' new_config
#'
#' @export
dmain_config <- list(
  r=10,
  lambda_1=1.0,
  lambda_2=1.0,
  n_iter=1000,
  converge=1e-6
)


#' @param X data matrix (p genes x n spots), input data
#' @param G graph matrix, input data
#' @param dmain_config hyper-parameters
#' @export
dmain <- function(X, G, dmain_config){

  ##check hyper-parameters
  config <- dmain_check(dmain_config)
  out <- fast(Xr=as.matrix(X), Gr=G, r=config$r, lambda_1r=config$lambda_1,
              lambda_2r=config$lambda_2, n_iter=config$n_iter, converge=config$converge)
  out
}
