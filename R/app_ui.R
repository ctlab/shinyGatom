# as action button, but doesn't restore state (?)
myActionButton <- function (inputId, label, icon = NULL, ...) {
    tags$button(id = inputId, type = "button", class = "btn btn-default action-button",
                list(icon, label), ...)
}

mySidebarPanel <- function(...) {
    div(class="sidebar col-sm-3", tags$form(class="well", ...))
}

myMainPanel <- function(...) {
    div(class="mainPanel col-sm-9", ...)
}

#' @import shinyCyJS
app_ui <- function(config_file) {
  conf <- config::get(file=config_file, use_parent = FALSE)

  function(request) {
    workPanel <- tagList(
        fixedRow(
            mySidebarPanel(
                actionButton(
                    "resetInput", label="Reset all",
                    onclick=paste(
                        'resetFileInput("geneDE"); resetFileInput("metDE")',
                        '$("#loadExampleGeneDE").prop("checked", false).trigger("change")',
                        '$("#loadExampleMetDE").prop("checked", false).trigger("change")',
                        '$("#loadExampleLipidDE").prop("checked", false).trigger("change")',
                        '$("#preprocess").click(); $("#find").click()',
                        sep=";")
                ),
                conditionalPanel("!input.loadExampleLipidDE",
                                 checkboxInput("loadExampleGeneDE",
                                               label="Example DE for genes",
                                               value=FALSE),
                                 checkboxInput("loadExampleMetDE",
                                               label="Example DE for metabolites",
                                               value=FALSE)),
                conditionalPanel("!input.loadExampleGeneDE * !input.loadExampleMetDE",
                                 checkboxInput("loadExampleLipidDE",
                                               label="Example DE for lipids",
                                               value=FALSE)),
                conditionalPanel(
                    "input.loadExampleGeneDE || input.loadExampleMetDE",
                    p("Organism: Mouse")
                ),
                conditionalPanel("!input.loadExampleLipidDE",
                                 selectInput("network",
                                             label="Select type of network",
                                             choices=c("KEGG network"="kegg",
                                                       "Rhea network"="rhea",
                                                       "Rhea lipid subnetwork"="lipidomic"),
                                             selected="kegg")
                                 ),
                conditionalPanel("input.loadExampleLipidDE",
                                 p("Network: Lipidomic Rhea subnetwork")),
                conditionalPanel(
                    "(!input.loadExampleGeneDE && !input.loadExampleMetDE) && !input.loadExampleLipidDE",
                    selectInput("organism",
                                label="Select an organism",
                                choices=c("Mouse"="mmu",
                                          "Human"="hsa",
                                          "Arabidopsis"="ath",
                                          "Yeast"="sce"),
                                selected="mmu"
                    ),
                    fileInput("geneDE", "File with DE for genes"),
                    fileInput("metDE", "File with DE for metabolites")
                ),
                conditionalPanel("(!input.loadExampleLipidDE) && (input.network != 'lipidomic')",
                                 selectInput("nodesAs",
                                             label="Interpret nodes as",
                                             choices=c("metabolites"="metabolites",
                                                       "atoms"="atoms"),
                                             selected="atoms")),
                conditionalPanel("(input.loadExampleLipidDE) || (input.network == 'lipidomic')",
                                 p("Interpret nodes as: metabolites")),
                conditionalPanel(
                    "false",
                    checkboxInput("autoFindModule",
                                  label="Set FDRs & run step 2 automatically",
                                  value=TRUE),
                    actionButton("preprocess",
                                 label="Step 1: Make network")
                ),
                myActionButton(
                    "runAll",
                    label="Make network and run step 2",
                    onclick='$("#autoFindModule")[0].checked=true; $("#autoFindModule").trigger("change"); $("#preprocess").trigger("click")',
                    disabled=""),
                div("or", style="text-align: center"),
                myActionButton(
                    "runStep1",
                    label="Step 1: Make network",
                    onclick='$("#autoFindModule")[0].checked=false; $("#autoFindModule").trigger("change"); $("#preprocess").trigger("click")',
                    disabled=""
                    )
            ),
            myMainPanel(
                div(class="DEBlock",
                    h3("Differential expression for genes"),
                    uiOutput("geneDESummary"),
                    # uiOutput("geneDEColumns"),
                    uiOutput("geneDETable"),
                    uiOutput("geneDENotMapped"),
                    uiOutput("geneDENotMappedTable")
                    ),
                div(class="DEBlock",
                    h3("Differential expression for metabolites"),
                    uiOutput("metDESummary"),
                    uiOutput("metDETable"),
                    uiOutput("metDENotMapped"),
                    uiOutput("metDENotMappedTable")
                    ),
                div(id="network-panel", class="bottom-buffer",
                    h3("Network summary"),
                    uiOutput("networkSummary"),
                    conditionalPanel(
                        "network.available",
                        downloadButton("downloadNetwork", "Download XGMML"),
                        # textOutput("test"),
                        div(id="bum-plots",
                            conditionalPanel("!input.nullkgene && network.hasGenes",
                                             class="col-sm-6",
                                             h4("BUM distribution for genes"),
                                             plotOutput("genesBUMPlot", width = "400px", height = "400px")
                            ),
                            conditionalPanel("!input.nullkmet && network.hasMets",
                                             class="col-sm-6",
                                             h4("BUM distribution for metabolites"),
                                             plotOutput("metsBUMPlot", width = "400px", height = "400px")
                            )
                        )
                    )
                ),
            )
        ),
        div(id="module-panel", class="row top-buffer",
            div(class="sidebar col-sm-3",
                wellPanel(
                    conditionalPanel("network.hasGenes",
                                     checkboxInput("nullkgene",
                                                   label="Don't use genes for scoring",
                                                   value=FALSE),
                                     numericInput("kgene",
                                                  label=HTML("k.gene"),
                                                  max=100, min=0, value=25, step=1)),

                    conditionalPanel("network.hasMets",
                                     checkboxInput("nullkmet",
                                                   label="Don't use metabolites for scoring",
                                                   value=FALSE),
                                     numericInput("kmet",
                                                  label=HTML("k.met"),
                                                  max=100, min=0, value=25, step=1)
                    ),

                    checkboxInput(
                        "solveToOptimality",
                        label="Try to solve to optimality",
                        value=FALSE),
                    conditionalPanel("network.available",
                                     uiOutput("solverString")
                    ),
                    myActionButton("find", "Step 2: Find module", disabled=""),
                    conditionalPanel("network.available",
                                     checkboxInput("addHighlyExpressedEdges",
                                                   label="Add highly expressed genes",
                                                   value=FALSE),
                                     selectInput("metaboliteActions",
                                                 label="Actions with atoms",
                                                 c("Connect atoms inside metabolite"="connectAtomsInsideMetabolite",
                                                   "Collapse atoms into metabolites"="collapseAtomsIntoMetabolites",
                                                   "None" = "NoAction"),
                                                 selected="NoAction")
                    )
                ),
                conditionalPanel("module.available",
                                 h3("Module summary"),
                                 uiOutput("moduleSummary"),
                                 myActionButton(
                                     "downloadSvg",
                                     label="SVG",
                                     onclick='var svgContent = cy.svg({scale: 1, full: true});\
				var blob = new Blob([svgContent], {type:"image/svg+xml;charset=utf-8"});\
				saveAs(blob, "module.svg");',
                                     icon=shiny::icon("download")),
                                 downloadButton("downloadModule", "XGMML")
                ),
                div(id="legend",
                    p(),
                    div(img(src="www/img/log2FCscale.svg"))
                ),
                p(),
                downloadButton("downloadVizMap", "Cytoscape VizMap")
            ),
            myMainPanel(
                h3("Module"),
                ShinyCyJSOutput(outputId = 'module', width = "100%", height = "100vh")
            )
        )
    )

    helpPanel <- fixedRow(
        mainPanel(id="helpPanel",
                  includeMarkdown(system.file("help.markdown", package="ShinyGATOM"))))

    aboutPanel <- fixedRow(
        mainPanel(includeMarkdown(system.file("about.markdown", package="ShinyGATOM"))))

    fluidPage(
        fluidPage(
            tags$head(
                # tags$script(src = "www/cytoscape-panzoom.js"),
                # tags$link(rel = "stylesheet", type = "text/css", href = "www/cytoscape.js-panzoom.css"),
                tags$script(src="www/d3.v3.min.js"),
                #tags$script(src="www/d3.v3.js"),
                tags$script(src="www/gam.js"),
                tags$link(rel="stylesheet", href="www/gam.css"),
                includeScript(system.file("ga.js", package="ShinyGATOM")),
                tags$title("Shiny GATOM")
            ),

            includeHTML(system.file("misc.xhtml", package="ShinyGATOM")),
            div(id="updateEsParameters", class="js-output"),

            titlePanel("Shiny GATOM: integrated analysis of genes and metabolites"),
            div(id="initialized", class="js-output"),

            fixedRow(
                column(12,
                       tabsetPanel(
                           tabPanel("Work", workPanel),
                           tabPanel("Help", helpPanel),
                           tabPanel("About", aboutPanel)
                       )
                ))
        )
    )
  }
}
