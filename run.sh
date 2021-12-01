#!/bin/sh
exec R -e "options(shiny.autoload.r=FALSE); shiny::runApp('.', launch.browser=F, port=8081)"
