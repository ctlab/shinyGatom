#!/bin/sh
exec R -e "options(shiny.autoload.r=FALSE,shiny.error=traceback,shiny.trace=TRUE,shiny.fullstacktrace=TRUE); shiny::runApp('.', launch.browser=F, port=8081)"
