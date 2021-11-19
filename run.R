library(devtools)
library(logging)
# load_all("/home/masha/Research/RCode/RShiny/shinyCyJS_Space")
# devtools::install()

# addHandler(writeToFile, file="/dev/stderr", level='DEBUG')

options(shiny.error=traceback)
options(shiny.trace=TRUE)
options(shiny.fullstacktrace=TRUE)
# shiny::devmode(TRUE)

"%o%" <- pryr::compose

Sys.setenv("R_CONFIG_ACTIVE"="local")
Sys.setenv("GATOM_LOCAL_DATA"="/home/masha/Research/RCode/RShiny/data")
Sys.setenv("CPLEX_HOME"="/home/masha/cplex")


## Running app
load_all()
ShinyGATOM()


## Testing
load_all()
# debug(ShinyGATOM)
shinytest::recordTest()
