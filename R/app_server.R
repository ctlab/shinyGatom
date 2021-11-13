#' @import mwcsr
#' @import data.table
app_server <- function(config_file) {

    conf <- config::get(file=config_file, use_parent = FALSE)

    function(input, output, session) {

        network_type <- list(
            "kegg"="network.kegg",
            "rhea"="network.rhea"
        )



        v.solver <- virgo_solver(cplex_dir=Sys.getenv("CPLEX_HOME"), penalty=0.01, timelimit=240)
        attr(v.solver, "description") <- "Virgo Solver (time limit = 4m)"


        v2.solver <- virgo_solver(cplex_dir=NULL,  penalty=0.01, timelimit=30)
        #v2.solver <- virgo_solver(cplex_dir="$CPLEX_HOME", penalty=0.01, timelimit=30)
        attr(v2.solver, "description") <- "Virgo Solver (time limit = 30s)"


        # Path to example data
        #example.gene.de.path <- "http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"
        #example.met.de.path <- "http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv"
        #example.gene.de.path <- "data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"
        #example.met.de.path <- "data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"

        # Loading example data
        #example.gene.de <- force(as.data.table(read.table(text=getURL(example.gene.de.path), stringsAsFactors=FALSE, header=1)))

        example.gene.de <- force(fread(conf$example.gene.de.path))
        attr(example.gene.de, "name") <- basename(conf$example.gene.de.path)

        #fread
        #example.met.de <- force(as.data.table(read.table(text=getURL(example.met.de.path), stringsAsFactors=FALSE, header=1)))

        example.met.de <- force(fread(conf$example.met.de.path))
        attr(example.met.de, "name") <- basename(conf$example.met.de.path)

        longProcessStart <- function() {
            session$sendCustomMessage(type='showWaitMessage', list(value=T))
        }

        longProcessStop <- function() {
            session$sendCustomMessage(type='showWaitMessage', list(value=F))
        }

        output$initialized <- renderJs('$("#initializing").hide()')

        queryValues <- reactiveValues()

        observe({
            query <- parseQueryString(session$clientData$url_search)
            if ("organism" %in% names(query)) {
                updateSelectInput(session, "network", selected=query$organism)
                #values$queryOrganism <- query$organism
            }

            if ("geneDE_key" %in% names(query)) {
                geneDE_key <- gsub("[^a-z0-9]", "", query$geneDE_key)
                loginfo("found key: %s", geneDE_key)
                if (geneDE_key == "") {
                    geneDE_key <- NULL
                }
                queryValues$geneDE_key <- geneDE_key
            }
        })

        loadExample <- reactive({
            input$loadExampleGeneDE || input$loadExampleMetDE
        })

        getNetwork <- reactive({
            #net <- if (loadExample()) {
            #    "kegg.mouse.network"
            #} else {
            #    networks[[input$network]]
            #}

            #lazyLoad(network_type[[input$network_type]])

            if (input$network_type == "kegg"){
                path.to.network <- conf$path.to.kegg.network
            } else {
                path.to.network <- conf$path.to.rhea.network
            }

            res <- lazyReadRDS(network_type[[input$network_type]],
                               path = path.to.network)
            res
        })


        geneDEInput <- reactive({
            if (loadExample()) {
                if (input$loadExampleGeneDE) {
                    return(example.gene.de)
                }
                return(NULL)
            }

            if (!is.null(input$geneDE) && !is(input$geneDE, "data.frame")) {
                # Value was reset
                return(NULL)
            }

            if (!is.null(input$geneDE)) {
                loginfo("GeneDE file:")
                loginfo(capture.output(str(input$geneDE)))
                loginfo("reading gene.de: %s", input$geneDE$name)
                loginfo("      from file: %s", input$geneDE$datapath)

                path <- input$geneDE$datapath
                deName <- input$geneDE$name
            } else if (!is.null(queryValues$geneDE_key)) {
                # User has not uploaded a file yet but provided value by key
                path <- file.path(uploaderPath, queryValues$geneDE_key)
                loginfo("GeneDE from key, path: %s", path)
                if (!file.exists(path)) {
                    return(NULL)
                }
                deName <- queryValues$geneDE_key
            } else {
                # User has not uploaded a file yet and we don't have a key
                return(NULL)
            }


            res <- read.table.smart.de.gene(path, idsList=gene.id.map)
            logdebug(capture.output(str(res)))
            if (!all(necessary.de.fields %in% names(res))) {
                loginfo("not all fields in DE file: %s", input$geneDE$datapath)
                if (grepl("xlsx?$", input$geneDE$name)) {
                    stop("We do not support excel files yet, please, use tab-separated files instead")
                } else{
                    stop(paste0("Genomic differential expression data should contain at least these fields: ",
                                paste(necessary.de.fields, collapse=", ")))
                }
            }
            attr(res, "name") <- deName
            res
        })

        geneIdsType <- reactive({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }


            org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                                conf$path.to.org.Mm.eg.gatom.anno)
            res <- getGeneDEMeta(data, org.gatom.anno=org.Mm.eg.gatom.anno)$idType

            if (length(res) != 1) {
                stop("Can't determine type of IDs for genes")
            }
            res
        })


        output$geneDESummary <- renderUI({
            gene.de <- geneDEInput()
            ids.type <- geneIdsType()
            if (is.null(gene.de)) {
                return("There is no genomic data")
            }

            div(
                HTML(
                    vector2html(c(
                        "name" = attr(gene.de, "name"),
                        "length" = nrow(gene.de),
                        "ID type" = ids.type
                    ))),
                p("Top DE genes:"))
        })

        output$geneDETable <- renderTable({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            format(as.data.frame(head(data[order(pval)])), digits=3)
        })

        geneMapsCreate <- reactive({
            org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                                path = path.to.org.Mm.eg.gatom.anno)

            entrez2refseq <- data.table(na.omit(org.Mm.eg.gatom.anno$mapFrom$RefSeq))
            colnames(entrez2refseq) <- c("Entrez", "RefSeq")
            entrez2ensembl <- data.table(na.omit(org.Mm.eg.gatom.anno$mapFrom$Ensembl))
            colnames(entrez2ensembl) <- c("Entrez", "Ensembl")
            entrez2symbol <- data.table(na.omit(org.Mm.eg.gatom.anno$mapFrom$Symbol))
            colnames(entrez2symbol) <- c("Entrez", "Symbol")

            entrez2name <- org.Mm.eg.gatom.anno$genes
            colnames(entrez2name) <- c("Entrez", "gene_symbol")

            res <- merge(entrez2name, entrez2refseq, all.x=T, by="Entrez")
            res <- merge(res, entrez2ensembl, all.x=T, by="Entrez")
            res <- merge(res, entrez2symbol, all.x=T, by="Entrez")

            res <- na.omit(res)
            res <- data.table(res)
            res
        })

        notMappedGenes <- reactive({
            org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno", path = path.to.org.Mm.eg.gatom.anno)
            geneIT <- geneIdsType()

            if (is.null(geneIT)) {
                return(NULL)
            }

            #if (geneIT == network$gene.ids) {
            if (geneIT == org.Mm.eg.gatom.anno$baseId) {
                return(NULL)
            }
            gene.id.map <- geneMapsCreate()
            notMapped <- setdiff(geneDEInput()$ID, gene.id.map[[geneIT]])
            notMapped
        })

        output$geneDENotMapped <- renderUI({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedGenes()
            network <- getNetwork()

            div(
                p(sprintf("Not mapped to %s: %s", org.Mm.eg.gatom.anno$baseId, length(notMapped))),
                if (length(notMapped) > 0) {
                    p("Top unmapped genes:",
                      a("show",
                        id="geneDENotMappedTableShowBtn",
                        href="javascript:void()",
                        onclick='$("#geneDENotMappedTable").show();$("#geneDENotMappedTableShowBtn").hide();$("#geneDENotMappedTableHideBtn").show()'),
                      a("hide",
                        id="geneDENotMappedTableHideBtn",
                        href="javascript:void()",
                        onclick='$("#geneDENotMappedTable").hide();$("#geneDENotMappedTableShowBtn").show();$("#geneDENotMappedTableHideBtn").hide()',
                        style="display:none"),
                      tag("script", '$("#geneDENotMappedTable").hide()')
                    )
                } else NULL)
        })

        output$geneDENotMappedTable <- renderTable({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedGenes()
            if (length(notMapped) == 0) {
                return(NULL)
            }

            data <- data[order(pval)]
            data <- data[ID %in% notMapped]
            format(as.data.frame(head(data, n=20)), digits=3)
        })

        metDEInput <- reactive({
            if (loadExample()) {
                if (input$loadExampleMetDE) {
                    return(example.met.de)
                }
                return(NULL)
            }


            if (is.null(input$metDE)) {
                # User has not uploaded a file yet
                return(NULL)
            }
            if (!is(input$metDE, "data.frame")) {
                # Value was reset
                return(NULL)
            }

            loginfo("MetDE file:")
            loginfo(capture.output(str(input$metDE)))
            loginfo("reading met.de: %s", input$metDE$name)
            loginfo("     from file: %s", input$metDE$datapath)

            res <- read.table.smart.de.met(input$metDE$datapath)
            logdebug(capture.output(str(res)))
            if (!all(necessary.de.fields %in% names(res))) {
                loginfo("not all fields in DE file: %s", input$metDE$datapath)
                if (grepl("xlsx?$", input$metDE$name)) {
                    stop("We do not support excel files yet, please, use tab-separated files instead")
                } else {
                    stop(paste0("Metabolic differential expression data should contain at least these fields: ",
                                paste(necessary.de.fields, collapse=", ")))
                }
            }
            attr(res, "name") <- input$metDE$name
            res
        })

        metIdsType <- reactive({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }

            if (input$network_type == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db
            }

            res <- getMetDEMeta(data, met.db=met.db)$idType

            if (length(res) != 1) {
                stop("Can't determine type of IDs for metabolites")
            }
            res
        })


        metMapsCreate <- reactive({
            if (input$network_type == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db

                hmdb2kegg <- data.table(na.omit(met.db$mapFrom$HMDB))
                colnames(hmdb2kegg) <- c("HMDB", "KEGG")
                kegg2name <- met.db$metabolites[ , c("metabolite", "metabolite_name")]
                colnames(kegg2name) <- c("KEGG", "name")
                res <- merge(hmdb2kegg, kegg2name, all.x=T, by="KEGG")
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db

                hmdb2chebi <- data.table(na.omit(met.db$mapFrom$HMDB))
                colnames(hmdb2chebi) <- c("HMDB", "ChEBI")
                chebi2name <- met.db$metabolites[ , c("metabolite", "metabolite_name")]
                colnames(chebi2name) <- c("ChEBI", "name")
                res <- merge(hmdb2chebi, chebi2name, all.x=T, by="ChEBI")
            }
            res <- na.omit(res)
            res <- data.table(res)
            res
        })

        output$metDESummary <- renderUI({
            met.de <- metDEInput()
            ids.type <- metIdsType()
            if (is.null(met.de)) {
                return("There is no metabolic data")
            }

            div(
                HTML(
                    vector2html(c(
                        "name" = attr(met.de, "name"),
                        "length" = nrow(met.de),
                        "ID type" = ids.type
                    ))),
                p("Top DE metabolites:"))
        })

        output$metDETable <- renderTable({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            format(as.data.frame(head(data[order(pval)])), digits=3)
        })

        notMappedMets <- reactive({
            #network <- getNetwork()

            if (input$network_type == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db
            }

            metIT <- metIdsType()

            if (is.null(metIT)) {
                return(NULL)
            }

            if (metIT == met.db$baseId) {
                return(NULL)
            }
            met.id.map <- metMapsCreate()
            notMapped <- setdiff(metDEInput()$ID, met.id.map[[metIT]])
        })

        output$metDENotMapped <- renderUI({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedMets()
            network <- getNetwork()

            if (input$network_type == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db
            }

            div(
                p(sprintf("Not mapped to %s: %s", met.db$baseId, length(notMapped))),
                if (length(notMapped) > 0) {
                    p("Top unmapped metabolites:",
                      a("show",
                        id="metDENotMappedTableShowBtn",
                        href="javascript:void()",
                        onclick='$("#metDENotMappedTable").show();$("#metDENotMappedTableShowBtn").hide();$("#metDENotMappedTableHideBtn").show()'),
                      a("hide",
                        id="metDENotMappedTableHideBtn",
                        href="javascript:void()",
                        onclick='$("#metDENotMappedTable").hide();$("#metDENotMappedTableShowBtn").show();$("#metDENotMappedTableHideBtn").hide()',
                        style="display:none"),
                      tag("script", '$("#metDENotMappedTable").hide()')
                    )
                } else NULL)
        })

        output$metDENotMappedTable <- renderTable({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedMets()
            if (length(notMapped) == 0) {
                return(NULL)
            }

            data <- data[order(pval)]
            data <- data[ID %in% notMapped]
            format(as.data.frame(head(data, n=20)), digits=3)
        })

        output$updateEsParameters <- renderJs({
            selected <- if (!is.null(metDEInput())) "edges" else "nodes"
            return(sprintf("$('#reactionsAs')[0].selectize.setValue('%s')", selected))
        })

        experimentTag <- reactive({
            geneData <- geneDEInput()
            metData <- metDEInput()
            name <- if (!is.null(geneData)) attr(geneData, "name") else attr(metData, "name")
            tag <- name
            tag <- gsub("\\.gz$", "", tag)
            tag <- gsub("\\.([ct]sv|txt)$", "", tag)
            tag
        })
        #
        #     output$reactionsAsHolder <- renderUI({
        #         gene.de <- geneDEInput()
        #
        #         met.de <- metDEInput()
        #
        #         selected <- if (!is.null(met.de)) "edges" else "nodes"
        #
        #         selectInput("reactionsAs",
        #                     label="Interpret reactions as",
        #                       c("edges"="edges", "nodes"="nodes"),
        #                       selected=selected)
        #     })


        esInput <- reactive({
            input$preprocess
            loginfo("Preprocessing")
            #input$runAll
            network <- isolate(getNetwork())
            gene.de <- isolate(geneDEInput())
            gene.ids <- isolate(geneIdsType())

            met.de <- isolate(metDEInput())
            met.ids <- isolate(metIdsType())
            tag <- isolate(experimentTag())

            if (is.null(gene.de) && is.null(met.de)) {
                return(NULL)
            }

            longProcessStart()

            tryCatch({
                if (!is.null(gene.de)) {
                    gene.de <- gene.de[which(gene.de$pval < 1),]
                }

                if (!is.null(met.de)) {
                    met.de <- met.de[which(met.de$pval < 1),]
                }

                reactions.as.edges = isolate(input$reactionsAs) == "edges"

                if (input$network_type == "kegg"){
                    met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                               path = conf$path.to.met.kegg.db)
                    met.db <- met.kegg.db
                    met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID
                } else if (input$network_type == "rhea") {
                    met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                               path = conf$path.to.met.rhea.db)
                    met.db <- met.rhea.db
                    met.to.filter = NULL
                }

                topology = isolate(input$nodesAs)

                es <- makeMetabolicGraph(network=network,
                                         topology=topology,
                                         org.gatom.anno=org.Mm.eg.gatom.anno,
                                         gene.de=geneDEInput(),
                                         met.db=met.db,
                                         met.de=metDEInput(),
                                         met.to.filter=met.to.filter)

                # g$has.genes <-
                #es$tag <- tag
                attr(es, "tag") <- tag
                es
            }, finally=longProcessStop())
        })

        output$networkSummary <- reactive({
            g <- esInput()
            if (is.null(g)) {
                return("There is no built network")
            }

            vector2html(c(
                "number of nodes" = length(V(g)),
                "number of edges" = length(E(g))
            ))
        })


        kMet <- reactive({
            if (!input$nullkmet) {
                input$kmet
            } else {
                NULL
            }
        })

        kGene <- reactive({
            if (!input$nullkgene) {
                input$kgene
            } else {
                NULL
            }
        })

        gene.de <- isolate(geneDEInput())
        met.de <- isolate(metDEInput())

        if (!is.null(gene.de) && !is.null(met.de)) {
            k.gene <- isolate(kGene())
            k.met <- isolate(kMet())
        } else if (!is.null(gene.de)) {
            k.gene <- isolate(kGene())
            k.met <- NULL
        } else if (!is.null(met.de)) {
            k.gene <- NULL
            k.met <- isolate(kMet())
        } else {
            k.gene <- NULL
            k.met <- NULL
        }


        metBUM <- reactive({
            g <- esInput()
            k.met <- isolate(kMet())
            res <- fitMetsToBUM(g=g, k.met=k.met)
            res
        })


        genBUM <- reactive({
            g <- esInput()
            if (is.null(g)) {
                return(NULL)
            }
            k.gene <- isolate(kGene())
            res <- fitGenesToBUM(g=g, k.gene=k.gene)
            res
        })


        output$genesBUMPlot <- renderPlot({
            g <- esInput()
            if (is.null(g)) {
                return(NULL)
            }
            gen.bum <- isolate(genBUM())
            plot(gen.bum)
            #list(plot(gen.bum))
        })



        output$metsBUMPlot <- renderPlot({
            g <- esInput()
            if (is.null(g)) {
                return("There is no built network")
            }
            met.bum <- isolate(metBUM())
            plot(met.bum)
        })


        output$networkParameters <- reactive({
            es <- NULL
            tryCatch({
                es <- esInput()
            }, error=function(e) {})

            if (is.null(es)) {
                return("")
            }

            gene.de <- isolate(geneDEInput())
            met.de <- isolate(metDEInput())

            has.genes <- FALSE
            has.mets <- FALSE

            if (!is.null(gene.de)) {
                has.genes <- TRUE
            }
            if (!is.null(met.de)) {
                has.mets <- TRUE
            }

            res <- paste0(
                makeJsAssignments(
                    network.available = TRUE,
                    #network.hasReactionsAsNodes = !es$reactions.as.edges,
                    #network.hasReactionsAsNodes = FALSE,
                    #network.hasReactionsAsEdges = es$reactions.as.edges,
                    #network.hasReactionsAsEdges = TRUE,
                    #network.hasGenes = !is.null(es$fb.rxn),
                    #network.hasGenes = has.genes,
                    network.hasGenes = TRUE,
                    #network.hasMets = !is.null(es$fb.met),
                    #network.hasMets = has.mets#,
                    network.hasMets = TRUE
                )
            )


            if (isolate(input$autoFindModule)) {
                res <- paste0(res, generateFDRs(es))
                res <- paste0(res, '$("#find").trigger("click");')
            }

            res <- paste0(res, '$("#find").removeAttr("disabled").addClass("btn-default");')
            res <- paste0(res, '$("#resetFDRs").removeAttr("disabled").addClass("btn-default");')
            res
        })

        output$enableMakeNetwork <- renderJs({
            res <- ""
            canRun <- FALSE
            tryCatch({
                geneDE <- geneDEInput()
                gIT <- geneIdsType()
                metDE <- metDEInput()
                mIT <- metIdsType()

                canRun <- !is.null(geneDE) || !is.null(metDE)
            }, error=function(e) {
                # if anything happened, not running
            })

            if (canRun) {
                res <- paste0(res, '$("#runStep1").removeAttr("disabled").addClass("btn-default");')
                res <- paste0(res, '$("#runAll").removeAttr("disabled").addClass("btn-default");')
            } else {
                res <- paste0(res, '$("#runStep1").attr("disabled", "disabled");')
                res <- paste0(res, '$("#runAll").attr("disabled", "disabled");')
            }
            res
        })

        output$showModulePanel <- renderJs({
            if (!is.null(esInput())) { return("mp = $('#module-panel'); mp[0].scrollIntoView();")
            }
            # return("mp = $('#module-panel'); mp.hide();")
            return("")
        })

        getSolver <- reactive({
            es <- esInput()
            if (is.null(es)) {
                return(NULL)
            }

            if (input$solveToOptimality) {
                #h.solver
                v.solver
                #} else if (es$reactions.as.edges && !is.null(es$fb.rxn)) {
                #g.solver
            } else {
                v2.solver
                #h2.solver
            }
        })

        output$solverString <- reactive({
            es <- esInput()
            solver <- getSolver()
            if (!is.null(solver)) {
                sprintf("Solver: %s", attr(solver, "description"))
            } else {
                ""
            }
        })

        esScoredInput <- reactive({
            input$find
            #met.fdr <- isolate(metFDR())
            #gene.fdr <- isolate(geneFDR())
            #absent.met.score <- isolate(input$absentMetScore)
            #absent.rxn.score <- isolate(input$absentRxnScore)
            #absent.rxn.score <- 0

            es <- isolate(esInput())

            if (is.null(es)) {
                return(NULL)
            }

            longProcessStart()
            loginfo(paste0(attr(es, "tag"),".mp"#, # min p-value
                           #if (es$reactions.as.edges) ".re" else ".rn",
                           #if (TRUE) ".re" else ".rn",
                           #".mf=", format(log10(met.fdr), digist=2),
                           #".rf=", format(log10(gene.fdr), digist=2),
                           #".ams=", absent.met.score
                           #, ".ars=", absent.rxn.score
            ))


            gene.de <- isolate(geneDEInput())
            met.de <- isolate(metDEInput())

            if (!is.null(gene.de) && !is.null(met.de)) {
                k.gene <- isolate(kGene())
                k.met <- isolate(kMet())
            } else if (!is.null(gene.de)) {
                k.gene <- isolate(kGene())
                k.met <- NULL
            } else if (!is.null(met.de)) {
                k.gene <- NULL
                k.met <- isolate(kMet())
            } else {
                k.gene <- NULL
                k.met <- NULL
            }

            met.bum <- isolate(metBUM())
            gen.bum <- isolate(genBUM())

            res <- scoreGraph2(g=es, k.gene=k.gene, k.met=k.met,
                               metabolite.bum=met.bum, gene.bum=gen.bum
            )


            res$description.string <- paste0(attr(es, "tag")#,
                                             #if (es$reactions.as.edges) ".re" else ".rn",
                                             #if (TRUE) ".re" else ".rn",
                                             #".mf=", format(log10(met.fdr), digist=2),
                                             #".rf=", format(log10(gene.fdr), digist=2),
                                             #".ams=", absent.met.score
                                             #, ".ars=", absent.rxn.score
            )
            res
        })

        rawModuleInput <- reactive({
            input$find

            esScored <- isolate(esScoredInput())

            if (is.null(esScored)) {
                return(NULL)
            }

            longProcessStart()
            tryCatch({
                solver <- isolate(getSolver())
                #res <- induced.subgraph(esScored$subnet.scored, V(esScored$subnet.scored)[adj(adj(1))])

                #res <- findModule(esScored, solver=solver)

                inst <- normalize_sgmwcs_instance(esScored,
                                                  nodes.weight.column = "score",
                                                  edges.weight.column = "score",
                                                  nodes.group.by = "signal",
                                                  edges.group.by = "signal",
                                                  group.only.positive = TRUE)

                res <- solve_mwcsp(solver, inst)
                res <- res$graph

                if (is.null(res) || length(V(res)) == 0) {
                    stop("No module found")
                }
                res$description.string <- esScored$description.string
                res
            }, finally=longProcessStop())
        })

        moduleInput <- reactive({
            module <- rawModuleInput()
            if (is.null(module)) {
                return(NULL)
            }

            # for consistency
            module <- remove.vertex.attribute(module, "score")
            module <- remove.edge.attribute(module, "score")

            es <- isolate(esInput())

            if (input$addHighlyExpressedEdges) {
                module <- addHighlyExpressedEdges(module, es)
            }

            #mets.actions = isolate(input$metaboliteActions)

            if (input$metaboliteActions == "connectAtomsInsideMetabolite") {
                module <- connectAtomsInsideMetabolite(module)
            } else if (input$metaboliteActions == "collapseAtomsIntoMetabolites") {
                module <- collapseAtomsIntoMetabolites(module)
            }

            module$description.string <- rawModuleInput()$description.string

            module
        })

        output$moduleSummary <- reactive({
            module <- moduleInput()
            if (is.null(module)) {
                return("There is no module yet")
            }

            vector2html(c(
                "number of nodes" = length(V(module)),
                "number of edges" = length(E(module))
            ))
        })

        output$moduleParameters <- reactive({
            m <- NULL
            tryCatch({
                m <- moduleInput()
            }, error=function(e) {})
            makeJsAssignments(
                module.available = !is.null(m)
            )
        })


        #output$table <- renderTable({
        #    module <- moduleInput()
        #    if (is.null(module)) {
        #        return(NULL)
        #    }
        #    module_e <- as_data_frame(module, what="edges")
        #    module_v <- as_data_frame(module, what="vertices")
        #    #save(module, file="example_module2.Rda")
        #    #module_e
        #    module_v
        #})

        getColumnNames <- function(dataInput){
            return(
                dataName = colnames(dataInput)
            )
        }

        prepareForshinyCyJS <- reactive({
            module <- moduleInput()
            if (is.null(module)) {
                return(NULL)
            }
            vertex.table <- as_data_frame(module, what="vertices")
            edge.table <- as_data_frame(module, what="edges")
            rownames(vertex.table) <- NULL
            save(module, file="module.rda")

            vertex.tooltip <- as.data.frame(lapply(colnames(vertex.table), function(x){
                lapply(vertex.table[x], function(y){
                    paste0("<b>", x, ":</b> ", y)
                })
            }))

            vertex.table$tooltip <- paste(vertex.tooltip$name, vertex.tooltip$metabolite,
                                          vertex.tooltip$label, vertex.tooltip$url,
                                          vertex.tooltip$pval, vertex.tooltip$HMDB,
                                          vertex.tooltip$baseMean, vertex.tooltip$logPval,
                                          vertex.tooltip$signal, vertex.tooltip$signalRank,
                                          vertex.tooltip$index, sep="<br>")

            nodes <- getDotNodeStyleAttributes2(vertex.table)
            edges <- getDotEdgeStyleAttributes2(edge.table)

            nodes = buildElems(nodes, type = 'Node')
            edges = buildElems(edges, type = 'Edge')

            obj = shinyCyJS(c(nodes, edges))
            obj
        })

        output$module <- renderShinyCyJS(prepareForshinyCyJS())

        ########################################################################
        output$downloadNetwork <- downloadHandler(
            filename = reactive({ paste0("network.", tolower(esInput()$network$organism), ".xgmml") }),
            content = function(file) {
                saveModuleToXgmml(esInput()$subnet, file=file, name=tolower(esInput()$network$organism))
            })

        output$downloadModule <- downloadHandler(
            filename = reactive({ paste0(moduleInput()$description.string, ".xgmml") }),
            content = function(file) {
                saveModuleToXgmml(moduleInput(), file=file, moduleInput()$description.string)
            })

        dotFile <- reactive({
            m <- moduleInput()
            if (is.null(m)) {
                return(NULL)
            }
            longProcessStart()
            res <- tempfile(pattern="module", fileext=".dot")
            saveModuleToDot(m, file=res, name=m$description.string)
            #saveModuleToDot(m, file=res, name="whatevs")
            longProcessStop()
            res
        })

        svgFile <- reactive({
            df <- dotFile()
            if (is.null(df)) {
                return(NULL)
            }
            res <- paste0(df, ".svg")
            system2("neato", c("-Tsvg",
                               "-o", res,
                               df), stderr=NULL)
            res
        })

        pdfFile <- reactive({
            df <- dotFile()
            if (is.null(df)) {
                return(NULL)
            }
            res <- paste0(df, ".pdf")
            system2("neato", c("-Tpdf",
                               "-o", res,
                               df), stderr=NULL)
            res
        })

        output$downloadXlsx <- downloadHandler(
            filename = reactive({ paste0(moduleInput()$description.string, ".xlsx") }),
            content = function(file) {
                module <- moduleInput()
                es <- isolate(esScoredInput())
                stopifnot(require(xlsx))

                wb <- createWorkbook()

                vTable <- data.table(get.vertex.attributes(module))
                eTable <- data.table(get.edge.attributes(module, include.ends=T))
                #if (es$reactions.as.edges) {
                if (TRUE) {
                    metTable <- vTable
                    metTable[, nodeType := NULL]
                    rxnTable <- eTable
                } else {
                    metTable <- vTable[nodeType == "met",]
                    metTable[, nodeType := NULL]
                    rxnTable <- vTable[nodeType == "rxn",]
                    rxnTable[, nodeType := NULL]
                }

                metTable <- removeNAColumns(metTable)
                rxnTable <- removeNAColumns(rxnTable)

                addDataFrame(metTable, createSheet(wb, "metabolites"), row.names=F)
                addDataFrame(rxnTable, createSheet(wb, "reactions"), row.names=F)

                metInModule <- metTable$name
                geneInModule <- rxnTable$origin

                saveWorkbook(wb, file)
            })

        output$downloadPDF <- downloadHandler(
            filename = reactive({ paste0(moduleInput()$description.string, ".pdf") }),
            content = function(file) {
                file.copy(pdfFile(), file)
            })

        output$downloadDot <- downloadHandler(

            filename = reactive({ paste0(moduleInput()$description.string, ".dot") }),
            content = function(file) {
                file.copy(dotFile(), file)
            })

        output$downloadVizMap <- downloadHandler(
            filename = "GAM_VizMap.xml",
            content = function(file) {
                file.copy(
                    from=system.file("GAM_VizMap.xml", package="GAM"),
                    to=file)
            })

        output$GAMVersion <- renderUI({
            p(paste("GAM version:", sessionInfo()$otherPkgs$GAM$Revision))
        })

        output$geneDEExample <- renderUI({
            a("here", href=config::get("example.gene.de.path", config="default", file=config_file, use_parent = FALSE), target="_blank")
        })

        output$metDEExample <- renderUI({
            a("here", href=config::get("example.met.de.path", config="default", file=config_file, use_parent = FALSE), target="_blank")
        })


        output$downloadVizMapInHelp <- downloadHandler(
            filename = "GAM_VizMap.xml",
            content = function(file) {
                file.copy(
                    from=system.file("GAM_VizMap.xml", package="GAM"),
                    to=file)
            })
    }
}
