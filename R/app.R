#' starts Shiny GATOM application
#' @import shiny
#' @import data.table
#' @import gatom
#' @import igraph
#' @import "BioNet"
#' @import mwcsr
#' @import RCurl
#' @import parallel
#' @import pryr
#' @import fgsea
#' @import markdown
#' @export
shinyGatom <- function(config_file=system.file('config.yml', package = 'shinyGatom')) {
    addResourcePath('/www',
                    system.file('www', package = 'shinyGatom'))

    shinyApp(app_ui(config_file), app_server(config_file))
}
