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
                    conditionalPanel(
                        "(!input.loadExampleGeneDE && !input.loadExampleMetDE) && !input.loadExampleLipidDE",
                        selectInput("organism",
                                    label="Select an organism",
                                    choices=c("Mouse"="mmu",
                                              "Human"="hsa",
                                              "Arabidopsis"="ath",
                                              "Yeast"="sce"
                                    ),
                                    selected="mmu"
                        ),
                        fileInput("geneDE", "File with DE for genes"),
                        fileInput("metDE", "File with DE for metabolites")
                    ),
                    conditionalPanel("!input.loadExampleLipidDE",
                                     selectInput("network",
                                                 label="Network type",
                                                 choices=c("KEGG network"="kegg",
                                                           "Rhea network"="rhea",
                                                           "Rhea lipid subnetwork"="lipidomic"),
                                                 selected="kegg")
                    ),
                    conditionalPanel("input.loadExampleLipidDE",
                                     p("Network type: Lipidomic Rhea subnetwork")),
                    conditionalPanel("(!input.loadExampleLipidDE) && (input.network != 'lipidomic')",
                                     selectInput("nodesAs",
                                                 label="Network topology",
                                                 choices=c("metabolites"="metabolites",
                                                           "atoms"="atoms"),
                                                 selected="atoms")),
                    conditionalPanel("(input.loadExampleLipidDE) || (input.network == 'lipidomic')",
                                     p("Network topology: metabolites")),
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
                    )
                )
            ),
            div(id="scoring-panel", class="row top-buffer",
                div(class="sidebar col-sm-3",
                    wellPanel(
                        conditionalPanel("network.hasGenes",
                                         checkboxInput("nullkgene",
                                                       label="Don't use genes for scoring",
                                                       value=FALSE),
                                         conditionalPanel("!input.nullkgene",
                                                          selectInput("scoringParameterGenes",
                                                                      label="Scoring parameter for genes",
                                                                      choices=c("Number of positive genes"="kgenechosen",
                                                                                "P-value threshold"="thresholdgenechosen",
                                                                                "FDR"="fdrgenechosen"),
                                                                      selected="kgenechosen"),

                                                          conditionalPanel('input.scoringParameterGenes == "kgenechosen"',
                                                                           numericInput("kgene",
                                                                                        label=HTML("Number of positive genes"),
                                                                                        max=100, min=0, value=50, step=1)),
                                                          conditionalPanel('input.scoringParameterGenes == "thresholdgenechosen"',
                                                                           numericInput("thresholdgene",
                                                                                        label=HTML("log<sub>10</sub> p-value threshold for genes"),
                                                                                        max=0, min=-100, value=-6, step=0.5)),
                                                          conditionalPanel('input.scoringParameterGenes == "fdrgenechosen"',
                                                                           numericInput("fdrgene",
                                                                                        label=HTML("log<sub>10</sub> FDR for genes"),
                                                                                        max=0, min=-100, value=-6, step=0.5))
                                                          )
                        ),

                        conditionalPanel("network.hasMets",
                                         checkboxInput("nullkmet",
                                                       label="Don't use metabolites for scoring",
                                                       value=FALSE),
                                         conditionalPanel("!input.nullkmet",
                                                          selectInput("scoringParameterMets",
                                                                      label="Scoring parameter for metabolites",
                                                                      choices=c("Number of positive metabolites"="kmetchosen",
                                                                                "P-value threshold"="thresholdmetchosen",
                                                                                "FDR"="fdrmetchosen"),
                                                                      selected="kmetchosen"),

                                                          conditionalPanel('input.scoringParameterMets == "kmetchosen"',
                                                                           numericInput("kmet",
                                                                                        label=HTML("Number of positive metabolites"),
                                                                                        max=100, min=0, value=50, step=1)),
                                                          conditionalPanel('input.scoringParameterMets == "thresholdmetchosen"',
                                                                           numericInput("thresholdmet",
                                                                                        label=HTML("log<sub>10</sub> p-value threshold for metabolites"),
                                                                                        max=0, min=-100, value=-6, step=0.5)),
                                                          conditionalPanel('input.scoringParameterMets == "fdrmetchosen"',
                                                                           numericInput("fdrmet",
                                                                                        label=HTML("log<sub>10</sub> FDR for metabolites"),
                                                                                        max=0, min=-100, value=-6, step=0.5)))
                        ),

                        checkboxInput(
                            "solveToOptimality",
                            label="Try to solve to optimality",
                            value=FALSE),

                        conditionalPanel("network.available",
                                         uiOutput("solverString")
                        ),

                        myActionButton("find", "Step 2: Find module", disabled="")
                    ),
                ),
                myMainPanel(
                    h3("Network summary"),
                    uiOutput("networkSummary"),
                    conditionalPanel(
                        "network.available",
                        downloadButton("downloadNetwork", "Download XGMML"),
                        # # :todo: rework the UI and uncomment
                        h3("Scoring"),
                        div(id="bum-plots",
                            conditionalPanel("network.hasGenes",
                                             class="col-sm-6",
                                             h4("BUM distribution for genes"),
                                             conditionalPanel("input.nullkgene",
                                                              p("Genes are not used for scoring")),
                                             conditionalPanel("!input.nullkgene",
                                                              uiOutput("genesBUMSummary"),
                                                              plotOutput("genesBUMPlot", width = "400px", height = "400px"),
                                                              uiOutput("genesScoringSummary")
                                             )
                            ),
                            conditionalPanel("network.hasMets",
                                             class="col-sm-6",
                                             h4("BUM distribution for metabolites"),
                                             conditionalPanel("input.nullkmet",
                                                              p("Metabolites are not used for scoring")),
                                             conditionalPanel("!input.nullkmet",
                                                              uiOutput("metsBUMSummary"),
                                                              plotOutput("metsBUMPlot", width = "400px", height = "400px"),
                                                              uiOutput("metsScoringSummary"))
                            )
                        )
                    )
                )
            ),
            div(id="module-panel", class="row top-buffer",
                div(class="sidebar col-sm-3",
                    wellPanel(
                        h4("Postprocessing"),
                        checkboxInput("addHighlyExpressedEdges",
                                      label="Add highly expressed genes",
                                      value=FALSE),
                        conditionalPanel("(input.nodesAs == 'atoms') * (!input.loadExampleLipidDE) * (input.network != 'lipidomic')",
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
                                     downloadButton("downloadModule", "XGMML"),
                                     downloadButton("downloadXlsx", "XLSX")
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
                    conditionalPanel("module.available",
                                     myActionButton("findPathways", "Show main pathways")
                    ),
                    ShinyCyJSOutput(outputId = 'module', width = "100%", height = "100vh")
                )
            ),
            conditionalPanel("pathways.show",
                         div(id="pathway-panel", class="row top-buffer",
                             div(class="sidebar col-sm-3",
                                 h3(""),
                                 wellPanel(
                                     checkboxInput("collapsePathways",
                                                   label="Collapse pathways",
                                                   value=FALSE)
                                 )
                             ),
                             myMainPanel(
                                 h3("Pathway annotation"),
                                 uiOutput("pathwaysTable")
                             )
                         )
        )
    )

        helpPanel <- fixedRow(
            mainPanel(id="helpPanel",
                      includeMarkdown(system.file("help.markdown", package="shinyGatom"))))

        aboutPanel <- fixedRow(
            mainPanel(includeMarkdown(system.file("about.markdown", package="shinyGatom"))))

        fluidPage(
            fluidPage(
                tags$head(
                    # tags$script(src = "www/cytoscape-panzoom.js"),
                    # tags$link(rel="stylesheet", type="text/css", href="www/cytoscape.js-panzoom.css"),
                    tags$script(src="www/d3.v3.min.js"),
                    #tags$script(src="www/d3.v3.js"),
                    tags$script(src="www/gam.js"),
                    tags$link(rel="stylesheet", href="www/gam.css"),
                    includeScript(system.file("ga.js", package="shinyGatom")),
                    tags$title("Shiny GATOM")
                ),

                includeHTML(system.file("misc.xhtml", package="shinyGatom")),
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
                    )
                )
            )
        )
    }
}
