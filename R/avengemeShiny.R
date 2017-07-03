#' Launch Shiny app for avengeme
#'
#' Launches a Shiny graphical user interface for the avengeme package.
#'
#' Uses the \code{\link{estimateGeneticModel}} function to estimate all parameters that are not explicitly set in the graphical interface.
#' If all parameters are set, calls the \code{\link{polygenescore}} function instead.
#' Parameters for combined genetic and environmental risk scores are not yet implemented.
#'
#' @seealso \code{\link{estimatePolygenicModel}}, \code{\link{polygenescore}}
#'
#' @author Frank Dudbridge

#' @references
#' Palla L and Dudbridge F (2015) A fast method using polygenic scores to estimate the variance explained by genome-wide marker panels and the proportion of variants affecting a trait. Am J Hum Genet 97:250-259
#' @export
avengemeShiny <- function() {
  appDir <- system.file("shiny", "avengeme", package = "avengeme")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `avengeme`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
