#' Launch GUI for running AVENGEME functions.
#'
#' Requires shiny
#'
#' @import shiny
AVENGEME_gui <- function(){
	shiny::runApp(system.file("shiny", package = "AVENGEME"))
		      
}
	

