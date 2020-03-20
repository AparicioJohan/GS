

#' User interface for running GS
#'
#' @return shiny interface
#' @export
#'
#' @examples
#' # GS::GS_shiny()
#' @author Johan Aparicio, \email{j.aparicio@cgiar.org}
#'
GS_shiny <- function()
{
  if (!requireNamespace(package = "shiny"))
    message("Package 'shiny' is required to run this function")
  if (!requireNamespace(package = "plotly"))
    message("Package 'plotly' is required to run this function")
  shiny::shinyAppDir(system.file("examples/",
                                  package = "GS", mustWork = TRUE), options = list(launch.browser=T))
}
