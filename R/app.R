# load_all("/home/masha/Research/RCode/RShiny/shinyCyJS")


#' starts shiny GATOM application
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
ShinyGATOM <- function(config_file=system.file('config.yml', package = 'ShinyGATOM')) {
    addResourcePath('/www',
                    system.file('www', package = 'ShinyGATOM'))

    shinyApp(app_ui, app_server(config_file))
}
