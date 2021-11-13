# library(shiny)
library(devtools)
library(logging)
load_all("/home/masha/Research/RCode/RShiny/shinyCyJS_edge")
# devtools::install()

load_all()

addHandler(writeToFile, file="/dev/stderr", level='DEBUG')

options(shiny.error=browser)
options(shiny.trace=TRUE)
options(shiny.fullstacktrace=TRUE)
shiny::devmode(TRUE)

"%o%" <- pryr::compose

ShinyGATOM()
debug(ShinyGATOM)
