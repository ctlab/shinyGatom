#' @import mwcsr
#' @import data.table
#' @import logging
app_server <- function(config_file) {
    
    conf <- config::get(file=config_file, use_parent = FALSE)
    
    function(input, output, session) {
        
        networks <- list(
            "kegg" = "network.kegg",
            "rhea" = "network.rhea"
        )
        
        networkPaths <- list(
            "kegg" = conf$path.to.kegg.network,
            "rhea" = conf$path.to.rhea.network
        )
        
        annotations <- list(
            "mmu" = "org.Mm.eg.gatom.anno",
            "hsa" = "org.Hs.eg.gatom.anno"
            # "ath" = "ath",
            # "sce" = "sce"
        )
        
        annotationPaths <- list(
            "mmu" = conf$path.to.org.Mm.eg.gatom.anno,
            "hsa" = conf$path.to.org.Hs.eg.gatom.anno
            # "ath"="ath",
            # "sce"="sce"
        )
        
        v.solver <- virgo_solver(cplex_dir=Sys.getenv("CPLEX_HOME"), penalty=0.01, timelimit=240)
        attr(v.solver, "description") <- "Virgo Solver (time limit = 4m)"
        
        v2.solver <- virgo_solver(cplex_dir=NULL, timelimit=30)
        # v2.solver <- virgo_solver(cplex_dir=Sys.getenv("CPLEX_HOME"), penalty=0.01, timelimit=30)
        attr(v2.solver, "description") <- "Virgo Solver (time limit = 30s)"
        
        
        # Loading example data
        example.gene.de <- force(fread(conf$example.gene.de.path))
        attr(example.gene.de, "name") <- basename(conf$example.gene.de.path)
        
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
            res <- lazyReadRDS(networks[[input$network]],
                               path = networkPaths[[input$network]])
            res
        })
        
        getAnnotation <- reactive({
            if (loadExample()){
                res <- lazyReadRDS(annotations[["mmu"]],
                                   path=annotationPaths[["mmu"]])
            } else {
                res <- lazyReadRDS(annotations[[input$organism]],
                                   path=annotationPaths[[input$organism]])
            }
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
            
            org.gatom.anno <- getAnnotation()
            
            res <- getGeneDEMeta(data, org.gatom.anno=org.gatom.anno)$idType
            
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
            org.gatom.anno <- getAnnotation()
            entrez2refseq <- data.table(na.omit(org.gatom.anno$mapFrom$RefSeq))
            colnames(entrez2refseq) <- c("Entrez", "RefSeq")
            entrez2ensembl <- data.table(na.omit(org.gatom.anno$mapFrom$Ensembl))
            colnames(entrez2ensembl) <- c("Entrez", "Ensembl")
            entrez2symbol <- data.table(na.omit(org.gatom.anno$mapFrom$Symbol))
            colnames(entrez2symbol) <- c("Entrez", "Symbol")
            
            entrez2name <- org.gatom.anno$genes
            colnames(entrez2name) <- c("Entrez", "gene_symbol")
            
            res <- merge(entrez2name, entrez2refseq, all.x=T, by="Entrez")
            res <- merge(res, entrez2ensembl, all.x=T, by="Entrez")
            res <- merge(res, entrez2symbol, all.x=T, by="Entrez")
            
            res <- na.omit(res)
            res <- data.table(res)
            res
        })
        
        notMappedGenes <- reactive({
            org.gatom.anno <- getAnnotation()
            geneIT <- geneIdsType()
            
            if (is.null(geneIT)) {
                return(NULL)
            }
            
            if (geneIT == org.gatom.anno$baseId) {
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
            org.gatom.anno <- getAnnotation()
            
            div(
                p(sprintf("Not mapped to %s: %s", 
                          org.gatom.anno$baseId, 
                          length(notMapped))),
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
            
            if (input$network == "kegg"){
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
            ids.type <- metIdsType()
            
            if (input$network == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db
                
                data2kegg <- data.table(na.omit(met.db$mapFrom[[ids.type]]))
                colnames(data2kegg) <- c(ids.type, "KEGG")
                kegg2name <- met.db$metabolites[ , c("metabolite", "metabolite_name")]
                colnames(kegg2name) <- c("KEGG", "name")
                res <- merge(data2kegg, kegg2name, all.x=T, by="KEGG")
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db
                
                data2chebi <- data.table(na.omit(met.db$mapFrom[[ids.type]]))
                colnames(data2chebi) <- c(ids.type, "ChEBI")
                chebi2name <- met.db$metabolites[ , c("metabolite", "metabolite_name")]
                colnames(chebi2name) <- c("ChEBI", "name")
                res <- merge(data2chebi, chebi2name, all.x=T, by="ChEBI")
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
            
            if (input$network == "kegg"){
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
            notMapped
        })
        
        output$metDENotMapped <- renderUI({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedMets()
            network <- getNetwork()
            
            if (input$network == "kegg"){
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
        
        experimentTag <- reactive({
            geneData <- geneDEInput()
            metData <- metDEInput()
            name <- if (!is.null(geneData)) attr(geneData, "name") else attr(metData, "name")
            tag <- name
            tag <- gsub("\\.gz$", "", tag)
            tag <- gsub("\\.([ct]sv|txt)$", "", tag)
            tag
        })
        
        
        gInput <- reactive({
            input$preprocess
            loginfo("Preprocessing")
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
                    gene.de <- gene.de[which(gene.de$pval < 1), ]
                }
                
                if (!is.null(met.de)) {
                    met.de <- met.de[which(met.de$pval < 1), ]
                }
                
                if (input$network == "kegg"){
                    met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                               path = conf$path.to.met.kegg.db)
                    met.db <- met.kegg.db
                } else if (input$network == "rhea") {
                    met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                               path = conf$path.to.met.rhea.db)
                    met.db <- met.rhea.db
                }
                
                topology = isolate(input$nodesAs)
                org.gatom.anno <- getAnnotation()
                
                g <- makeMetabolicGraph(network=network,
                                         topology=topology,
                                         org.gatom.anno=org.gatom.anno,
                                         gene.de=geneDEInput(),
                                         met.db=met.db,
                                         met.de=metDEInput(),
                                         met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
                
                attr(g, "tag") <- tag
                g
            }, finally=longProcessStop())
        })
        
        output$networkSummary <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return("There is no built network")
            }
            
            vector2html(c(
                "number of nodes" = length(V(g)),
                "number of edges" = length(E(g))
            ))
        })
        
        
        
        kMet <- reactive({
            met.de <- isolate(metDEInput())
            
            if (is.null(met.de)) {
                return(NULL)
            } else if (input$nullkmet) {
                return(NULL)
            }
            input$kmet
        })
        
        kGene <- reactive({
            gene.de <- isolate(geneDEInput())
            
            if (is.null(gene.de)) {
                return(NULL)
            } else if (input$nullkgene) {
                return(NULL)
            }
            input$kgene
        })
        
        
        metBUM <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            k.met <- isolate(kMet())
            if (is.null(k.met)) {
                return(NULL)
            }
            res <- fitMetsToBUM(g=g, k.met=k.met)
            res
        })
        
        
        genBUM <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            k.gene <- isolate(kGene())
            if (is.null(k.gene)) {
                return(NULL)
            }
            res <- fitGenesToBUM(g=g, k.gene=k.gene)
            res
        })
        
        
        output$genesBUMPlot <- renderPlot({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            gen.bum <- isolate(genBUM())
            if (is.null(gen.bum)) {
                return(NULL)
            }
            
            plot(gen.bum)
        })
        
        
        output$metsBUMPlot <- renderPlot({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            met.bum <- isolate(metBUM())
            if (is.null(met.bum)) {
                return(NULL)
            }
            
            plot(met.bum)
        })
        
        
        output$networkParameters <- reactive({
            g <- NULL
            tryCatch({
                g <- gInput()
            }, error=function(e) {})
            
            if (is.null(g)) {
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
                    network.hasGenes = has.genes,
                    network.hasMets = has.mets
                )
            )
            
            
            if (isolate(input$autoFindModule)) {
                res <- paste0(res, '$("#find").trigger("click");')
            }
            
            res <- paste0(res, '$("#find").removeAttr("disabled").addClass("btn-default");')
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
            if (!is.null(gInput())) { return("mp = $('#network-panel'); mp[0].scrollIntoView();")
            }
            return("")
        })
        
        getSolver <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            
            if (input$solveToOptimality) {
                v.solver
            } else {
                v2.solver
            }
        })
        
        output$solverString <- reactive({
            g <- gInput()
            solver <- getSolver()
            if (!is.null(solver)) {
                sprintf("Solver: %s", attr(solver, "description"))
            } else {
                ""
            }
        })
        
        gScoredInput <- reactive({
            input$find
            g <- isolate(gInput())
            
            if (is.null(g)) {
                return(NULL)
            }
            
            longProcessStart()
            loginfo(paste0(attr(g, "tag"),".mp"
            ))
            
            
            gene.de <- isolate(geneDEInput())
            met.de <- isolate(metDEInput())
            
            k.gene <- isolate(kGene())
            k.met <- isolate(kMet())
            
            met.bum <- isolate(metBUM())
            gen.bum <- isolate(genBUM())
            
            res <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                                   metabolite.bum=met.bum, gene.bum=gen.bum
            )
            
            
            res$description.string <- paste0(attr(g, "tag")#,
                                             #if (g$reactions.as.edges) ".re" else ".rn",
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
            
            gScored <- isolate(gScoredInput())
            
            if (is.null(gScored)) {
                return(NULL)
            }
            
            longProcessStart()
            tryCatch({
                solver <- isolate(getSolver())
                
                inst <- gScored
                
                res <- solve_mwcsp(solver, inst)
                res <- res$graph
                
                if (is.null(res) || length(V(res)) == 0) {
                    stop("No module found")
                }
                res$description.string <- gScored$description.string
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
            
            g <- isolate(gInput())
            
            if (input$addHighlyExpressedEdges) {
                module <- addHighlyExpressedEdges(module, g)
            }
            
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
        
        output$module <- renderShinyCyJS(prepareForShinyCyJS(moduleInput()))

        output$downloadNetwork <- downloadHandler(
            filename = reactive({ paste0("network.", tolower(gInput()$network$organism), ".xgmml") }),
            content = function(file) {
                saveModuleToXgmml(gInput()$subnet, file=file, name=tolower(gInput()$network$organism))
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
