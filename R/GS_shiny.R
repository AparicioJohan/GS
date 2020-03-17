

#' user interface for running GS
#'
#' @return
#' @export
#'
#' @examples
GS_shiny <- function()
{
  if (!requireNamespace(package = "shiny"))
    message("Package 'shiny' is required to run this function")
  if (!requireNamespace(package = "plotly"))
    message("Package 'plotly' is required to run this function")
  shiny::shinyAppDir(system.file("examples/",
                                  package = "GS", mustWork = TRUE))
}
