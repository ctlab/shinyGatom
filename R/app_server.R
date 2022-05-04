#' @import mwcsr
#' @import data.table
#' @import logging
#' @import openxlsx
#' @import shinyjs
app_server <- function(config_file) {

    conf <- config::get(file=config_file, use_parent = FALSE)

    function(input, output, session) {

        networks <- list(
            "kegg" = "network.kegg",
            "rhea" = "network.rhea",
            "lipidomic" = "network.rhea.lipids"
        )

        networkPaths <- list(
            "kegg" = conf$path.to.kegg.network,
            "rhea" = conf$path.to.rhea.network,
            "lipidomic" = conf$path.to.lipid.network
        )

        annotations <- list(
            "mmu" = "org.Mm.eg.gatom.anno",
            "hsa" = "org.Hs.eg.gatom.anno",
            "ath" = "org.At.tair.gatom.anno",
            "sce" = "org.Sc.sgd.gatom.anno"
        )

        annotationPaths <- list(
            "mmu" = conf$path.to.org.Mm.eg.gatom.anno,
            "hsa" = conf$path.to.org.Hs.eg.gatom.anno,
            "ath" = conf$path.to.org.At.tair.gatom.anno,
            "sce" = conf$path.to.org.Sc.sgd.gatom.anno
        )

        annotationsRhea <- list(
            "mmu" = "gene2reaction.rhea.mmu.eg",
            "hsa" = "gene2reaction.rhea.hsa.eg",
            "ath" = "gene2reaction.rhea.ath.eg",
            "sce" = "gene2reaction.rhea.sce.eg"
        )

        annotationRheaPaths <- list(
            "mmu" = conf$path.to.gene2reaction.rhea.mmu.eg,
            "hsa" = conf$path.to.gene2reaction.rhea.hsa.eg,
            "ath" = conf$path.to.gene2reaction.rhea.ath.eg,
            "sce" = conf$path.to.gene2reaction.rhea.sce.eg
        )

        v.solver <- virgo_solver(cplex_dir=Sys.getenv("CPLEX_HOME"), penalty=0.01, timelimit=240, threads=4)
        attr(v.solver, "description") <- "Virgo Solver (time limit = 4m)"

        # v2.solver <- virgo_solver(cplex_dir=NULL, timelimit=30)
        v2.solver <- virgo_solver(cplex_dir=Sys.getenv("CPLEX_HOME"), penalty=0.01, timelimit=30, threads=4)
        attr(v2.solver, "description") <- "Virgo Solver (time limit = 30s)"


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
            if (input$loadExampleLipidDE) {
                res <- lazyReadRDS(networks[["lipidomic"]],
                                   path = networkPaths[["lipidomic"]])
            } else {
                res <- lazyReadRDS(networks[[input$network]],
                                   path = networkPaths[[input$network]])
            }
            res
        })

        getAnnotation <- reactive({
            if (loadExample()){
                res <- lazyReadRDS(annotations[["mmu"]],
                                   path=annotationPaths[["mmu"]])
            } else if (input$loadExampleLipidDE) {
                res <- lazyReadRDS(annotations[["mmu"]],
                                   path=annotationPaths[["mmu"]])
            } else {
                res <- lazyReadRDS(annotations[[input$organism]],
                                   path=annotationPaths[[input$organism]])
            }
            res
        })



        geneDEInputRaw <- reactive({
            if (loadExample()) {
                if (input$loadExampleGeneDE) {
                    example.gene.de <- force(fread(conf$example.gene.de.path))
                    attr(example.gene.de, "name") <- basename(conf$example.gene.de.path)
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

            if (grepl("xlsx?$", input$geneDE$name)){
                res <- read.xlsx(path)
            } else {
                res <- fread(path, colClasses="character")
            }
            attr(res, "name") <- deName

            res
        })


        geneDEMeta <- reactive({
            gene.de.raw <- geneDEInputRaw()

            if (is.null(gene.de.raw)) {
                return(NULL)
            }

            org.gatom.anno <- getAnnotation()
            gene.de.meta <- getGeneDEMeta(gene.de.raw, org.gatom.anno)
            gene.de.meta
        })

        geneDEInput <- reactive({
            gene.de.raw <- geneDEInputRaw()

            if (is.null(gene.de.raw)) {
                return(NULL)
            }

            org.gatom.anno <- getAnnotation()
            gene.de.meta <- geneDEMeta()
            loginfo(capture.output(str(gene.de.meta)))

            res <- prepareDE(gene.de.raw, gene.de.meta)
            # res[, signalRank := NULL] # todo

            logdebug(capture.output(str(res)))
            if (!all(necessary.gene.de.fields %in% names(res))) {
                loginfo("not all fields in DE file: %s", input$geneDE$datapath)
                stop(paste0("Genomic differential expression data should contain at least these fields: ",
                            paste(necessary.gene.de.fields, collapse=", ")))
            }

            attr(res, "name") <- attr(gene.de.raw, "name")

            res
        })


        geneIdsType <- reactive({
            gene.de.meta <- geneDEMeta()

            if (is.null(gene.de.meta)) {
                return(NULL)
            }

            res <- gene.de.meta$idType


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

            if ("signalRank" %in% colnames(data)){
                data[ , signalRank := NULL]
            }
            format(as.data.frame(head(data[order(pval)])), digits=3)
        })

        notMappedGenes <- reactive({
            org.gatom.anno <- getAnnotation()
            geneIT <- geneIdsType()
            gene.de <- geneDEInput()

            if (is.null(geneIT)) {
                return(NULL)
            }

            gene.de.meta <- geneDEMeta()
            data <- gene.de
            data$initialID <- data$ID

            split <- " */// *"

            if (length(split) != 0) {
                newIds <- gatom:::splitIds(data[["ID"]])
                if (any(lengths(newIds) != 1)) {
                    data.new <- data[rep(seq_len(nrow(data)),
                                         lengths(newIds))]
                    data.new[, `:=`(c("ID"), unlist(newIds))]
                    data <- data.new
                }
            }

            data[, c("ID") := sub(pattern="\\.\\d*$", replacement="", x=data[["ID"]])]

            if (geneIT == org.gatom.anno$baseId) {
                notMappedIDs <- setdiff(data$ID, org.gatom.anno$genes$gene)
            } else {
                notMappedIDs <- setdiff(data$ID, org.gatom.anno$mapFrom[[geneIT]][[geneIT]])
            }

            data.new <- data[!(data$ID %in% notMappedIDs)]
            data.new <- data.new[ , -"ID"]
            data.new <- data.new[!duplicated(data.new), ]

            data <- data[ , -"ID"]
            notMappedDT <- data[!(data$initialID %in% data.new$initialID)]
            notMapped <- notMappedDT$initialID

            notMapped
        })

        output$geneDENotMapped <- renderUI({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedGenes()
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
            not.mapped.data <- data[ID %in% notMapped]
            format(as.data.frame(head(not.mapped.data, n=20)), digits=3)
        })
        outputOptions(output, "geneDENotMappedTable", suspendWhenHidden = FALSE)
        
        
        geneDENotMappedTableToDownload <- reactive({
            data <- geneDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedGenes()
            if (length(notMapped) == 0) {
                return(NULL)
            }
            
            data <- data[order(pval)]
            not.mapped.data <- data[ID %in% notMapped]
            not.mapped.data
        })
        
        geneDETag <- reactive({
            geneData <- geneDEInput()
            tag <- attr(geneData, "name")
            tag <- gsub("\\.gz$", "", tag)
            tag <- gsub("\\.([ct]sv|txt|xlsx)$", "", tag)
            tag
        })
        
        output$inputDataParameters <- reactive({
            has.genes <- FALSE
            has.mets <- FALSE
            
            gene.de <- geneDEInput()
            unmapped.genes <- notMappedGenes()
            if (!is.null(unmapped.genes) && length(unmapped.genes) != 0) {
                has.genes <- TRUE
            }
            
            met.de <- metDEInput()
            unmapped.mets <- notMappedMets()
            if (!is.null(unmapped.mets) && length(unmapped.mets) != 0) {
                has.mets <- TRUE
            }
            
            res <- paste0(
                makeJsAssignments(
                    inputdata.hasGenes = has.genes,
                    inputdata.hasMets = has.mets
                )
            )
            res
        })
        
        output$downloadGeneDENotMappedTable <- downloadHandler(
            filename = reactive({ paste0(geneDETag(), ".unmapped.genes.csv")}),
            content = function(file) {
                write.csv(geneDENotMappedTableToDownload(), file, row.names = FALSE)
            }
        )

        metDEInputRaw <- reactive({
            if (loadExample()) {
                if (input$loadExampleMetDE) {
                    example.met.de <- force(fread(conf$example.met.de.path, colClasses="character"))
                    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
                    return(example.met.de)
                }
                return(NULL)
            }

            if (input$loadExampleLipidDE){
                example.lip.de <- force(fread(conf$example.lip.de.path, colClasses="character"))
                attr(example.lip.de, "name") <- basename(conf$example.lip.de.path)
                return(example.lip.de)
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

            if (grepl("xlsx?$", input$metDE$name)){
                res <- read.xlsx(input$metDE$datapath)
            } else {
                res <- fread(input$metDE$datapath, colClasses="character")
            }
            attr(res, "name") <- input$metDE$name

            res
        })
        
        getMetDb <- reactive({
            met.de.raw <- metDEInputRaw()
            if (is.null(met.de.raw)) {
                return(NULL)
            }
            
            if ((input$loadExampleLipidDE) || (input$network == "lipidomic")) {
                met.lipid.db <- lazyReadRDS(name = "met.lipid.db",
                                            path = conf$path.to.met.lipid.db)
                met.db <- met.lipid.db
            } else if (input$network == "kegg"){
                met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                           path = conf$path.to.met.kegg.db)
                met.db <- met.kegg.db
            } else {
                met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                           path = conf$path.to.met.rhea.db)
                met.db <- met.rhea.db
            }
            
            met.db
        })
        
        metDEMeta <- reactive({
            met.de.raw <- metDEInputRaw()
            if (is.null(met.de.raw)) {
                return(NULL)
            }
            
            met.db <- getMetDb()
            
            met.de.meta <- getMetDEMeta(met.de.raw, met.db)
            met.de.meta
        })

        metDEInput <- reactive({
            met.de.raw <- metDEInputRaw()
            if (is.null(met.de.raw)) {
                return(NULL)
            }
            
            met.db <- getMetDb()

            met.de.meta <- metDEMeta()

            res <- prepareDE(met.de.raw, met.de.meta)

            logdebug(capture.output(str(res)))
            if (!all(necessary.met.de.fields %in% names(res))) {
                loginfo("not all fields in DE file: %s", input$metDE$datapath)
                stop(paste0("Metabolic differential expression data should contain at least these fields: ",
                            paste(necessary.met.de.fields, collapse=", ")))
            }

            attr(res, "name") <- attr(met.de.raw, "name")
            res

        })


        metIdsType <- reactive({
            met.de.meta <- metDEMeta()
            if (is.null(met.de.meta)) {
                return(NULL)
            }

            res <- met.de.meta$idType

            if (length(res) != 1) {
                stop("Can't determine type of IDs for metabolites")
            }
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

            if ("signalRank" %in% colnames(data)){
                data[ , signalRank := NULL]
            }
            format(as.data.frame(head(data[order(pval)])), digits=3)
        })


        notMappedMets <- reactive({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            
            metIT <- metIdsType()
            if (is.null(metIT)) {
                return(NULL)
            }
            
            data$initialID <- data$ID

            split <- " */// *"

            if (length(split) != 0) {
                newIds <- gatom:::splitIds(data[["ID"]])
                if (any(lengths(newIds) != 1)) {
                    data.new <- data[rep(seq_len(nrow(data)),
                                         lengths(newIds))]
                    data.new[, `:=`(c("ID"), unlist(newIds))]
                    data <- data.new
                }
            }
            
            met.db <- getMetDb()
            
            if (metIT == met.db$baseId) {
                notMappedIDs <- setdiff(data$ID, met.db$metabolites$metabolite)
            } else {
                notMappedIDs <- setdiff(data$ID, met.db$mapFrom[[metIT]][[metIT]])
            }

            data.new <- data[!(data$ID %in% notMappedIDs)]
            data.new <- data.new[ , -"ID"]

            data.new <- data.new[!duplicated(data.new), ]

            data <- data[ , -"ID"]
            notMappedDT <- data[!(data$initialID %in% data.new$initialID)]
            notMapped <- notMappedDT$initialID

            notMapped
        })


        output$metDENotMapped <- renderUI({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedMets()
            met.db <- getMetDb()

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
        outputOptions(output, "metDENotMappedTable", suspendWhenHidden = FALSE)

        
        metDENotMappedTableToDownload <- reactive({
            data <- metDEInput()
            if (is.null(data)) {
                return(NULL)
            }
            notMapped <- notMappedMets()
            if (length(notMapped) == 0) {
                return(NULL)
            }
            
            data <- data[order(pval)]
            not.mapped.data <- data[ID %in% notMapped]
            not.mapped.data
        })
        
        metDETag <- reactive({
            metData <- metDEInput()
            tag <- attr(metData, "name")
            tag <- gsub("\\.gz$", "", tag)
            tag <- gsub("\\.([ct]sv|txt|xlsx)$", "", tag)
            tag
        })
        
        output$downloadMetDENotMappedTable <- downloadHandler(
            filename = reactive({ paste0(metDETag(), ".unmapped.metabolites.csv")}),
            content = function(file) {
                write.csv(metDENotMappedTableToDownload(), file, row.names = FALSE)
            }
        )
        
        
        experimentTag <- reactive({
            geneData <- geneDEInput()
            metData <- metDEInput()
            name <- if (!is.null(geneData)) attr(geneData, "name") else attr(metData, "name")
            tag <- name
            tag <- gsub("\\.gz$", "", tag)
            tag <- gsub("\\.([ct]sv|txt|xlsx)$", "", tag)
            tag
        })
        
        
        typeIsSpecies <- reactive({
            met.de <- metDEInput()
            if (is.null(met.de)) {
                return(FALSE)
            }
            
            metIT <- isolate(metIdsType())
            loginfo(metIT)
            
            if (metIT == "Species") {
                flag <- TRUE
            } else {
                flag <- FALSE
            }
            return(flag)
        })
        
        originalNamesFlag <- reactive({
            met.de <- metDEInput()
            if (is.null(met.de)) {
                return(FALSE)
            }
            speciesFlag <- typeIsSpecies()
            loginfo(speciesFlag)
            if (speciesFlag) {
                if (input$useOriginalLipidNames) {
                    return(TRUE)
                }
            }
            
            return(FALSE)
        })
        
        output$lipidSpeciesParameters <- reactive({
            flag <- typeIsSpecies()
            
            makeJsAssignments(
                mapping.fromSpecies = flag
            )
        })
        
        gInput <- reactive({
            input$preprocess
            loginfo("Preprocessing")

            gene.de <- isolate(geneDEInput())
            met.de <- isolate(metDEInput())

            if (is.null(gene.de) && is.null(met.de)) {
                return(NULL)
            }

            network <- isolate(getNetwork())
            gene.ids <- isolate(geneIdsType())
            met.ids <- isolate(metIdsType())
            tag <- isolate(experimentTag())


            longProcessStart()

            tryCatch({

                topology <- isolate(input$nodesAs)
                org.gatom.anno <- isolate(getAnnotation())
                keepReactionsWithoutEnzymes <- FALSE
                gene2reaction.extra <- NULL

                if ((isolate(input$network) == "lipidomic") || isolate(input$loadExampleLipidDE)) {
                    met.lipid.db <- lazyReadRDS(name = "met.lipid.db",
                                                path = conf$path.to.met.lipid.db)
                    met.db <- met.lipid.db

                    if (isolate(input$loadExampleLipidDE)) {
                        gene2reaction.extra <- (fread(annotationRheaPaths[["mmu"]], colClasses="character"))[gene != "-"]
                    } else {
                        gene2reaction.extra <- (fread(annotationRheaPaths[[input$organism]], colClasses="character"))[gene != "-"]
                    }
                    
                    if (met.ids == "Species") {
                        met.de$SpecialSpeciesLabelColumn <- met.de$ID
                    }
                    
                    topology <- "metabolites"
                    keepReactionsWithoutEnzymes <- TRUE
                } else if (isolate(input$network) == "kegg"){
                    met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                               path = conf$path.to.met.kegg.db)
                    met.db <- met.kegg.db
                } else {
                    met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                               path = conf$path.to.met.rhea.db)
                    met.db <- met.rhea.db
                }

                g <- makeMetabolicGraph(network=network,
                                        topology=topology,
                                        org.gatom.anno=org.gatom.anno,
                                        gene.de=gene.de,
                                        met.db=met.db,
                                        met.de=met.de,
                                        met.to.filter=fread(system.file("mets2mask.lst",
                                                                        package="gatom"))$ID,
                                        keepReactionsWithoutEnzymes = keepReactionsWithoutEnzymes,
                                        gene2reaction.extra = gene2reaction.extra)

                attr(g, "tag") <- tag
                g$organism <- isolate(input$organism)
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


        kGene <- reactive({
            input$preprocess
            gene.de <- isolate(geneDEInput())

            if (is.null(gene.de)) {
                return(NULL)
            }
            if (input$nullkgene) {
                return(NULL)
            }
            input$kgene
        })

        thresholdGene <- reactive({
            input$preprocess
            gene.de <- isolate(geneDEInput())

            if (is.null(gene.de)) {
                return(NULL)
            } else if (input$nullkgene) {
                return(NULL)
            }
            10^input$thresholdgene
        })

        fdrGene <- reactive({
            input$preprocess
            gene.de <- isolate(geneDEInput())

            if (is.null(gene.de)) {
                return(NULL)
            } else if (input$nullkgene) {
                return(NULL)
            }
            10^input$fdrgene
        })


        kMet <- reactive({
            input$preprocess
            met.de <- isolate(metDEInput())

            if (is.null(met.de)) {
                return(NULL)
            } else if (input$nullkmet) {
                return(NULL)
            }
            input$kmet
        })

        thresholdMet <- reactive({
            input$preprocess
            met.de <- isolate(metDEInput())

            if (is.null(met.de)) {
                return(NULL)
            } else if (input$nullkmet) {
                return(NULL)
            }
            10^input$thresholdmet
        })

        fdrMet <- reactive({
            input$preprocess
            met.de <- isolate(metDEInput())

            if (is.null(met.de)) {
                return(NULL)
            } else if (input$nullkmet) {
                return(NULL)
            }
            10^input$fdrmet
        })

        genBUM <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            gene.de <- isolate(geneDEInput())
            if (is.null(gene.de)) {
                return(NULL)
            }

            if (is.null(input$kgene) & is.null(input$thresholdgene) & is.null(input$fdrgene)) {
                return(NULL)
            }

            edge.table <- data.table(as_data_frame(g, what="edges"))

            res <- fitToBUM(table=edge.table)
            res
        })


        metBUM <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            met.de <- isolate(metDEInput())
            if (is.null(met.de)) {
                return(NULL)
            }

            if (is.null(input$kmet) & is.null(input$fdrmet) & is.null(input$thresholdmet)) {
                return(NULL)
            }

            vertex.table <- data.table(as_data_frame(g, what="vertices"))

            res <- fitToBUM(table=vertex.table)
            res
        })


        genesScoring <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            edge.table <- data.table(as_data_frame(g, what="edges"))

            gene.bum <- genBUM()
            if (is.null(gene.bum)) {
                return(NULL)
            }

            if (input$scoringParameterGenes == 'kgenechosen') {
                k.gene <- kGene()
                scores <- getScoringParameters(table=edge.table,
                                               scoring.parameter="k",
                                               scor.par.value=k.gene,
                                               bum=gene.bum)

            } else if (input$scoringParameterGenes == 'thresholdgenechosen') {
                threshold.gene <- thresholdGene()
                scores <- getScoringParameters(table=edge.table,
                                               scoring.parameter="threshold",
                                               scor.par.value=threshold.gene,
                                               bum=gene.bum)
            } else {
                fdr.gene <- fdrGene()
                scores <- getScoringParameters(table=edge.table,
                                               scoring.parameter="fdr",
                                               scor.par.value=fdr.gene,
                                               bum=gene.bum)
            }

            return(scores)
        })


        output$genesBUMSummary <- renderUI({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            edge.table <- data.table(as_data_frame(g, what="edges"))

            gene.bum <- isolate(genBUM())
            if (is.null(gene.bum)) {
                return(NULL)
            }

            pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

            div(
                p(sprintf("Gene BU alpha: %.3g", gene.bum$a)),
                p(sprintf("Amount of gene signals: %s", length(pvalsToFit)))
            )
        })

        output$genesBUMPlot <- renderPlot({

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            edge.table <- data.table(as_data_frame(g, what="edges"))

            gen.bum <- genBUM()
            if (is.null(gen.bum)) {
                return(NULL)
            }

            hist(gen.bum)
        })

        output$genesScoringSummary <- renderUI({
            req(input$kgene, input$thresholdgene, input$fdrgene)

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            gene.bum <- isolate(genBUM())
            if (is.null(gene.bum)) {
                return(NULL)
            }

            scores <- genesScoring()
            if (is.null(scores)) {
                return(NULL)
            }

            if (gene.bum$a > 0.5) {
                div(
                    p(sprintf("Edge scores have been assigned to 0 due to an inappropriate p-value distribution"))
                )
            } else {
                div(
                    p(sprintf("Number of positive genes: %s", round(scores$k, 0))),
                    p(sprintf("Gene p-value threshold: %.1g", scores$threshold)),
                    p(sprintf("FDR for genes: %.1g", scores$fdr))
                )
            }
        })


        metsScoring <- reactive({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            vertex.table <- data.table(as_data_frame(g, what="vertices"))

            met.bum <- isolate(metBUM())
            if (is.null(met.bum)) {
                return(NULL)
            }

            if (input$scoringParameterMets == 'kmetchosen') {
                k.met <- kMet()
                scores <- getScoringParameters(table=vertex.table,
                                               scoring.parameter="k",
                                               scor.par.value=k.met,
                                               bum=met.bum)

            } else if (input$scoringParameterMets == 'thresholdmetchosen') {
                threshold.met <- thresholdMet()
                scores <- getScoringParameters(table=vertex.table,
                                               scoring.parameter="threshold",
                                               scor.par.value=threshold.met,
                                               bum=met.bum)
            } else {
                fdr.met <- fdrMet()
                scores <- getScoringParameters(table=vertex.table,
                                               scoring.parameter="fdr",
                                               scor.par.value=fdr.met,
                                               bum=met.bum)
            }

            return(scores)
        })

        output$metsBUMSummary <- renderUI({
            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            vertex.table <- data.table(as_data_frame(g, what="vertices"))

            met.bum <- isolate(metBUM())
            if (is.null(met.bum)) {
                return(NULL)
            }

            pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

            div(
                p(sprintf("Metabolite BU alpha: %.3g", met.bum$a)),
                p(sprintf("Amount of metabolite signals: %s", length(pvalsToFit)))
            )
        })

        output$metsBUMPlot <- renderPlot({

            req(input$kmet | input$thresholdmet | input$fdrmet)

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }
            met.bum <- metBUM()
            if (is.null(met.bum)) {
                return(NULL)
            }

            hist(met.bum)
        })

        output$metsScoringSummary <- renderUI({
            req(input$kmet, input$thresholdmet, input$fdrmet)

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            met.bum <- isolate(metBUM())
            if (is.null(met.bum)) {
                return(NULL)
            }

            scores <- metsScoring()
            if (is.null(scores)) {
                return(NULL)
            }

            if (met.bum$a > 0.5) {
                div(
                    p(sprintf("Node scores have been assigned to 0 due to an inappropriate p-value distribution"))
                )
            } else {
                div(
                    p(sprintf("Number of positive metabolites: %s", round(scores$k, 0))),
                    p(sprintf("Metabolite p-value threshold: %.1g", scores$threshold)),
                    p(sprintf("FDR for metabolites: %.1g", scores$fdr))
                )
            }
        })



        observeEvent(input$preprocess, {

            output$networkParameters <- reactive({
                g <- NULL
                tryCatch({
                    g <- isolate(gInput())
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

        output$showScoringPanel <- renderJs({
            if (!is.null(gInput())) { return("$('#scoring-panel')[0].scrollIntoView()")
            }
            return("")
        })

        output$showModulePanel <- renderJs({
            if (!is.null(moduleInput())) { return("$('#module-panel')[0].scrollIntoView()")
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

            gene.bum <- genBUM()
            metabolite.bum <- metBUM()

            gen.scores <- genesScoring()
            met.scores <- metsScoring()

            res <- scoreGraphShiny(g=g, k.gene=gen.scores$k, k.met=met.scores$k,
                                   metabolite.bum=metabolite.bum, gene.bum=gene.bum,
                                   vertex.threshold=met.scores$threshold,
                                   edge.threshold=gen.scores$threshold
            )


            res$description.string <- paste0(attr(g, "tag")
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

        output$module <- renderShinyCyJS(prepareForShinyCyJS(moduleInput(), 
                                                             orig_names=originalNamesFlag()))

        output$downloadNetwork <- downloadHandler(
            filename = reactive({ paste0("network.", tolower(gInput()$organism), ".xgmml")}),
            content = function(file) {
                saveModuleToXgmml(gInput(), file=file, name=tolower(gInput()$organism))
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

        output$downloadXlsx <- downloadHandler(
            filename = reactive({ paste0(moduleInput()$description.string, ".xlsx") }),
            content = function(file) {
                module <- moduleInput()
                g <- isolate(gScoredInput())

                stopifnot(require(openxlsx))

                wb <- createWorkbook()

                metTable <- data.table(as_data_frame(module, what="vertices"))
                rxnTable <- data.table(as_data_frame(module, what="edges"))

                addWorksheet(wb, "metabolites")
                writeData(wb, "metabolites", metTable, row.names=F)

                addWorksheet(wb, "reactions")
                writeData(wb, "reactions", rxnTable, row.names=F)

                metInModule <- metTable$name
                geneInModule <- rxnTable$gene

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
                vizmap.url <- conf$vizmap.url
                if (startsWith(vizmap.url, "https://") || startsWith(vizmap.url, "http://")) {
                    download.file(vizmap.url, destfile=file)
                } else {
                    file.copy(
                        from=vizmap.url,
                        to=file)
                }
            })

        output$GatomVersion <- renderUI({
            p(paste("gatom version:", packageVersion("gatom")))
        })


        observeEvent(input$findPathways, {

            output$pathwaysParameters <- reactive({
                    makeJsAssignments(
                        pathways.show = TRUE
                    )
                })
        })

        output$showPathwaysPanel <- renderJs({
            if (input$findPathways) {
                return("setTimeout(() => $('#pathway-panel')[0].scrollIntoView(), 50)")
            }
            return("")
        })


        pathwaysInput <- reactive({
            m <- moduleInput()
            if (is.null(m)) {
                return(NULL)
            }

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            org.gatom.anno <- getAnnotation()

            foraRes <- fgsea::fora(pathways=org.gatom.anno$pathways,
                                   genes=E(m)$gene,
                                   universe=unique(E(g)$gene),
                                   minSize=5)
            res <- foraRes[padj < 0.05]
#
#             if (input$collapsePathways) {
#                 mainPathways <- fgsea::collapsePathwaysORA(
#                     res,
#                     pathways=org.gatom.anno$pathways,
#                     genes=E(m)$gene,
#                     universe=unique(E(g)$gene))
#
#                 res <- foraRes[pathway %in% mainPathways$mainPathways]
#             }

            res
        })


        output$pathwaysTable <- renderTable({
            m <- moduleInput()
            if (is.null(m)) {
                return(NULL)
            }

            g <- gInput()
            if (is.null(g)) {
                return(NULL)
            }

            res <- pathwaysInput()
            format(as.data.frame(res), digits=3)
        })

        observe({
            pathways <- pathwaysInput()

            pathwayChoices <- c("None", pathways$pathway)
            pathwayChoices <- as.character(pathwayChoices)
            names(pathwayChoices) <- pathwayChoices

            updateSelectInput(session, "selectPathway",
                              choices = pathwayChoices,
                              selected = "None")
        })

        pathwaySelection <- reactive({
            input$selectPathway
        })

        observeEvent(input$selectPathway, {
            selectedPath <- pathwaySelection()

            if (selectedPath != "None") {
                js$UnSelectPathway()

                m <- moduleInput()
                if (is.null(m)) {
                    return(NULL)
                }

                edge.table <- as_data_frame(m, what="edges")
                edge.table <- as.data.table(edge.table)
                edge.table <- edge.table[, gene:=as.character(gene)]
                edge.table <- as.data.frame(edge.table)

                pathways <- pathwaysInput()
                data <- pathways[pathways$pathway == selectedPath]
                overlapGenes <- unlist(data$overlapGenes)
                # ids <- row.names(edge.table[edge.table$gene %in% overlapGenes, ])
                # ids <- as.numeric(ids)

                genesToSelect <- edge.table[edge.table$gene %in% overlapGenes, ]
                genesToSelect <- genesToSelect$label

                Edges <- paste0('[label = "', genesToSelect, '"]', collapse = ",")
                js$SelectPathway(Edges)

            } else {
                js$UnSelectPathway()
            }
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
                vizmap.url <- conf$vizmap.url
                if (startsWith(vizmap.url, "https://") || startsWith(vizmap.url, "http://")) {
                    download.file(vizmap.url, destfile=file)
                } else {
                    file.copy(
                        from=vizmap.url,
                        to=file)
                }
            })
    }
}
