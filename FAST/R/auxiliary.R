# package FAST
# functions for argument checking


#' Funtion to check dmain arguments 
#'
#' @keywords internal
#' @noRd
#' @param config list with dmain arguments
#'
#' @return config object
dmain_check = function(config=dmain_config, ...) {
  ##these values should be integers
  config$r <- ceiling(config$r)
  config$n_iter <- ceiling(config$n_iter)
  if(config$r < 2) {
    stop(paste0("rank must be greater than 1", "\n"), call.=FALSE)
  }
  if(config$n_iter < 1) {
    stop(paste0("number of iterations must be greater than 0", "\n"), call.=FALSE)
  }
  config
}


