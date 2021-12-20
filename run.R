library(devtools)
library(logging)

options(shiny.error=traceback)
options(shiny.trace=TRUE)
options(shiny.fullstacktrace=TRUE)
# shiny::devmode(TRUE)

Sys.setenv("R_CONFIG_ACTIVE"="local")
Sys.setenv("GATOM_LOCAL_DATA"="/home/masha/Research/RCode/RShiny/data")
Sys.setenv("CPLEX_HOME"="/home/masha/cplex")


## Debugging
# debug(ShinyGATOM)


## Running app
load_all()
shinyGatom()


## Testing
# load_all()
# devtools::test()
# shinytest::recordTest()
