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
#' @import markdown
ShinyGATOM <- function() {
    addResourcePath('/www',
                    system.file('www', package = 'ShinyGATOM'))

    #source("./functions.R")
    #source("./config.R")

    #data("met.id.map")
    #data("kegg.human.network")
    #data("kegg.mouse.network")
    #data("kegg.arabidopsis.network")
    #data("kegg.yeast.network")

    #networks <- list(
    #    "mmu"="kegg.mouse.network",
    #    "hsa"="kegg.human.network",
    #    "ath"="kegg.arabidopsis.network",
    #    "sce"="kegg.yeast.network"
    #    )

    shinyApp(app_ui, app_server)
}
