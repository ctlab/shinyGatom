
getJsTooltip <- function(attr.values) {
    attr.strings <- list()
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.strings[[i]] <- sprintf("<b>%s:</b> %s", attr.names[i], attr.values[[i]])
    }
    names(attr.strings) <- attr.names
    apply(do.call("cbind", attr.strings), 1, paste0, collapse="<br>")
}


getJsNodeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(
        width=sapply(logPval, gatom:::getDotSize) * 60,
        height=sapply(logPval, gatom:::getDotSize) * 60,
        label=if (!is.null(attrs$label)) label else "",
        id=attrs$name,
        shape=if (!is.null(attrs$nodeType)) nodeShapeMap[nodeType] else "ellipse",
        fontSize= getFontSizeJs(logPval), # fontSize=sapply(logPval, gatom:::getDotSize) * 40,
        bgColor=if (!is.null(attrs$log2FC)) sapply(log2FC, gatom:::getDotColor) else "#7777ff",
        borderWidth=4,
        borderColor="#eee",
        labelColor="black",
        tooltip=getJsTooltip(attrs)
    ))
}

getFontSizeJs <- function(val) {
    val <- sapply(val, gatom:::getDotSize) * 40
    val <- replace(val, val < 15, 15)
    val
}

getJsEdgeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(
        source=attrs$from,
        target=attrs$to,
        label=if (!is.null(attrs$label)) label else if (!is.null(attrs$gene)) gene else "",
        lineStyle="solid",
        fontSize= getFontSizeJs(logPval), # sapply(logPval, gatom:::getDotSize) * 40,
        lineColor=if (!is.null(attrs$log2FC)) sapply(as.numeric(log2FC), gatom:::getDotColor) else "grey",
        tooltip=getJsTooltip(attrs)
    ))
}

# #' @import gatom
# scoreGraphShiny <- function(g, k.gene, k.met,
#                             vertex.threshold.min=0.1,
#                             edge.threshold.min=0.1,
#                             met.score.coef=1,
#                             metabolite.bum=NULL,
#                             gene.bum=NULL,
#                             show.warnings=TRUE,
#                             raw=FALSE) {
#     if (show.warnings) {
#         warnWrapper <- identity
#     } else {
#         warnWrapper <- suppressWarnings
#     }
# 
#     vertex.table <- data.table(as_data_frame(g, what="vertices"))
#     edge.table <- data.table(as_data_frame(g, what="edges"))
# 
#     if (!is.null(k.met)) {
#         pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
#         if(is.null(metabolite.bum)) {
#             warnWrapper(metabolite.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
#         }
# 
#         if (metabolite.bum$a > 0.5) {
#             V(g)$score <- 0
#             warning("Vertex scores have been assigned to 0 due to an inappropriate p-value distribution")
# 
#         } else {
#             vertex.threshold <- if (k.met > length(pvalsToFit)) 1 else {
#                 sort(pvalsToFit)[k.met]
#             }
# 
#             vertex.threshold <- min(vertex.threshold,
#                                     BioNet::fdrThreshold(vertex.threshold.min, metabolite.bum))
# 
#             met.fdr <- gatom:::.reversefdrThreshold(vertex.threshold, metabolite.bum)
# 
#             gatom:::.messagef("Metabolite p-value threshold: %f", vertex.threshold)
#             gatom:::.messagef("Metabolite BU alpha: %f", metabolite.bum$a)
#             gatom:::.messagef("FDR for metabolites: %f", met.fdr)
# 
#             V(g)$score <- with(vertex.table,
#                                (metabolite.bum$a - 1) *
#                                    (log(gatom:::.replaceNA(pval, 1)) - log(vertex.threshold)))
#             V(g)$score <- V(g)$score * met.score.coef
#         }
#     }
#     else {
#         V(g)$score <- 0
#         V(g)$signal <- ""
# 
#     }
# 
#     if (!is.null(k.gene)) {
#         pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
#         if(is.null(gene.bum)){
#             warnWrapper(gene.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
#         }
# 
#         if(gene.bum$a > 0.5) {
#             E(g)$score <- 0
#             warning("Edge scores have been assigned to 0 due to an inappropriate p-value distribution")
# 
#         } else {
#             edge.threshold <- if (k.gene > length(pvalsToFit)) 1 else {
#                 sort(pvalsToFit)[k.gene]
#             }
# 
#             edge.threshold <- min(edge.threshold,
#                                   BioNet::fdrThreshold(edge.threshold.min, gene.bum))
# 
#             gene.fdr <- gatom:::.reversefdrThreshold(edge.threshold, gene.bum)
# 
#             gatom:::.messagef("Gene p-value threshold: %f", edge.threshold)
#             gatom:::.messagef("Gene BU alpha: %f", gene.bum$a)
#             gatom:::.messagef("FDR for genes: %f", gene.fdr)
# 
#             E(g)$score <- with(edge.table,
#                                (gene.bum$a - 1) *
#                                    (log(gatom:::.replaceNA(pval, 1)) - log(edge.threshold)))
#         }
#     }
#     else {
#         E(g)$score <- 0
#         E(g)$signal <- ""
#     }
#     g
#     if (raw) {
#         return(g)
#     }
# 
#     res <- normalize_sgmwcs_instance(g,
#                                      nodes.weight.column = "score",
#                                      edges.weight.column = "score",
#                                      nodes.group.by = "signal",
#                                      edges.group.by = "signal",
#                                      group.only.positive = TRUE)
#     res
# }

#' @import gatom
getScoringParameters <- function(table, 
                                 scoring.parameter=c("fdr", "k", "threshold"), 
                                 scor.par.value,
                                 bum=NULL,
                                 threshold.min=0.1,
                                 show.warnings=TRUE,
                                 raw=FALSE){
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }
    
    if (is.null(scor.par.value)) {
        return(list(threshold = NULL,
                    k = NULL,
                    fdr = NULL,
                    signals = NULL))
    }
    
    pvalsToFit <- table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
    signals <- length(pvalsToFit)
    
    if(is.null(bum)){
        warnWrapper(bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
    }
    
    if (bum$a > 0.5) {
        return(list(threshold = NULL,
                    k = NULL,
                    fdr = NULL,
                    signals = signals))
    }
    
    if (scoring.parameter == "k") {
        k <- scor.par.value
        threshold <- if (k > length(pvalsToFit)) 1 else {
            sort(pvalsToFit)[k]
        }
        
        threshold <- min(threshold,
                         BioNet::fdrThreshold(threshold.min, bum))
    }
    
    if (scoring.parameter != "fdr"){
        if (scoring.parameter == "threshold") {
            threshold <- scor.par.value
        }
        fdr <- gatom:::.reversefdrThreshold(threshold, bum)
    } else {
        fdr <- scor.par.value
        threshold <- BioNet::fdrThreshold(fdr, bum)
    }
    
    if (scoring.parameter != "k"){
        if (threshold > max(pvalsToFit)) {
            k <- length(pvalsToFit) + 1
        } else {
            pvalsToFit <- sort(pvalsToFit)
            k <- length(pvalsToFit[pvalsToFit <= threshold])
        }
    }
    
    res <- list(threshold = threshold,
                k = k,
                fdr = fdr,
                signals = signals)
    
    return(res)
}



# #' @import gatom
# scoreNetworkFromK <- function(table, 
#                               k,
#                               bum=NULL,
#                               threshold.min=0.1,
#                               show.warnings=TRUE,
#                               raw=FALSE){
#     if (show.warnings) {
#         warnWrapper <- identity
#     } else {
#         warnWrapper <- suppressWarnings
#     }
#     
#     if (!is.null(k)) {
#         pvalsToFit <- table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
#         signals <- length(pvalsToFit)
#         
#         if(is.null(bum)){
#             warnWrapper(bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
#         }
#         
#         if (bum$a <= 0.5) {
#             threshold <- if (k > length(pvalsToFit)) 1 else {
#                 sort(pvalsToFit)[k]
#             }
#             
#             threshold <- min(threshold,
#                              BioNet::fdrThreshold(threshold.min, bum))
#             
#             fdr <- gatom:::.reversefdrThreshold(threshold, bum)
#             
#             res <- list(threshold = threshold,
#                         k = k,
#                         fdr = fdr,
#                         signals = signals)
#         } else {
#             res <- list(threshold = NULL,
#                         k = NULL,
#                         fdr = NULL,
#                         signals = signals)
#         }
#     } else {
#         res <- list(threshold = NULL,
#                     k = NULL,
#                     fdr = NULL,
#                     signals = NULL)
#     }
#     return(res)
# }


# #' @import gatom
# scoreNetworkFromThreshold <- function(table, 
#                                       threshold,
#                                       bum=NULL,
#                                       threshold.min=0.1,
#                                       show.warnings=TRUE,
#                                       raw=FALSE){
#     if (show.warnings) {
#         warnWrapper <- identity
#     } else {
#         warnWrapper <- suppressWarnings
#     }
#     
#     if (is.null(threshold)) {
#         return(list(threshold = NULL,
#                     k = NULL,
#                     fdr = NULL,
#                     signals = NULL))
#     }
#     
#     pvalsToFit <- table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
#     signals <- length(pvalsToFit)
#     
#     if(is.null(bum)){
#         warnWrapper(bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
#     }
#     
#     if (bum$a > 0.5) {
#         return(list(threshold = NULL,
#                     k = NULL,
#                     fdr = NULL,
#                     signals = signals))
#     }
#     
#     fdr <- gatom:::.reversefdrThreshold(threshold, bum)
#     
#     if (threshold > max(pvalsToFit)) {
#         k <- length(pvalsToFit) + 1
#     } else {
#         pvalsToFit <- sort(pvalsToFit)
#         k <- length(pvalsToFit[pvalsToFit <= threshold])
#     }
#     
#     res <- list(threshold = threshold,
#                 k = k,
#                 fdr = fdr,
#                 signals = signals)
#     
#     
#     return(res)
# }


# #' @import gatom
# scoreNetworkFromFDR <- function(table, 
#                                 fdr,
#                                 bum=NULL,
#                                 threshold.min=0.1,
#                                 show.warnings=TRUE,
#                                 raw=FALSE){
#     if (show.warnings) {
#         warnWrapper <- identity
#     } else {
#         warnWrapper <- suppressWarnings
#     }
#     
#     if (!is.null(fdr)) {
#         pvalsToFit <- table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
#         signals <- length(pvalsToFit)
#         
#         if(is.null(bum)){
#             warnWrapper(bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
#         }
#         
#         if (bum$a <= 0.5) {
#             # threshold <- if (k > length(pvalsToFit)) 1 else {
#             #     sort(pvalsToFit)[k]
#             # }
#             
#             threshold <- BioNet::fdrThreshold(fdr, bum)
#             
#             # fdr <- gatom:::.reversefdrThreshold(threshold, bum)
#             
#             if (threshold > max(pvalsToFit)) {
#                 k <- length(pvalsToFit) + 1
#             } else {
#                 pvalsToFit <- sort(pvalsToFit)
#                 k <- length(pvalsToFit[pvalsToFit <= threshold])
#             }
#             
#             res <- list(threshold = threshold,
#                         k = k,
#                         fdr = fdr,
#                         signals = signals)
#         } else {
#             res <- list(threshold = NULL,
#                         k = NULL,
#                         fdr = NULL,
#                         signals = signals)
#         }
#     } else {
#         res <- list(threshold = NULL,
#                     k = NULL,
#                     fdr = NULL,
#                     signals = NULL)
#     }
#     return(res)
# }

#' @import gatom
scoreGraphShiny <- function(g, k.gene, k.met,
                            vertex.threshold=NULL,
                            edge.threshold=NULL,
                            met.score.coef=1,
                            metabolite.bum=NULL,
                            gene.bum=NULL,
                            show.warnings=TRUE,
                            raw=FALSE) {
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }

    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    edge.table <- data.table(as_data_frame(g, what="edges"))

    if (!is.null(k.met)) {
        if (metabolite.bum$a > 0.5) {
            V(g)$score <- 0
        } else {
            V(g)$score <- with(vertex.table,
                               (metabolite.bum$a - 1) *
                                   (log(gatom:::.replaceNA(pval, 1)) - log(vertex.threshold)))
            V(g)$score <- V(g)$score * met.score.coef
        }
    } else {
        V(g)$score <- 0
        V(g)$signal <- ""
    }
    
    
    if (!is.null(k.gene)) {
        if(gene.bum$a > 0.5) {
            E(g)$score <- 0
        } else {
            E(g)$score <- with(edge.table,
                               (gene.bum$a - 1) *
                                   (log(gatom:::.replaceNA(pval, 1)) - log(edge.threshold)))
        }
    } else {
        E(g)$score <- 0
        E(g)$signal <- ""
    }
    g
    
    if (raw) {
        return(g)
    }
    # delete
    # save(g, file="../error_in_normalizing.rda")

    res <- normalize_sgmwcs_instance(g,
                                     nodes.weight.column = "score",
                                     edges.weight.column = "score",
                                     nodes.group.by = "signal",
                                     edges.group.by = "signal",
                                     group.only.positive = TRUE)
    res
}


############################################### NEW FUNCTIONS


prepareForShinyCyJS <- function(module, layout=list(name="cose", animate=FALSE), orig_names){
    if (is.null(module)) {
        return(NULL)
    } else {

        vertex.table <- as_data_frame(module, what="vertices")
        edge.table <- as_data_frame(module, what="edges")
        
        if (orig_names) {
            vertex.table$tmp <- vertex.table$label
            vertex.table$label <- vertex.table$SpecialSpeciesLabelColumn
            vertex.table$SpecialSpeciesLabelColumn <- vertex.table$tmp
            vertex.table$tmp <- NULL
            
            vertex.table$label[is.na(vertex.table$label)] <- vertex.table$Abbreviation_LipidMaps[is.na(vertex.table$label)]
            vertex.table$label[is.na(vertex.table$label)] <- vertex.table$Abbreviation_SwissLipids[is.na(vertex.table$label)]
            vertex.table$label[vertex.table$label == "-"] <- vertex.table$Abbreviation_SwissLipids[vertex.table$label == "-"]
            vertex.table$label[is.na(vertex.table$label)] <- vertex.table$SpecialSpeciesLabelColumn[is.na(vertex.table$label)]
            vertex.table$label[vertex.table$label == "-"] <- vertex.table$SpecialSpeciesLabelColumn[vertex.table$label == "-"]
        }
        
        nodes <- getJsNodeStyleAttributes(vertex.table)
        edges <- getJsEdgeStyleAttributes(edge.table)

        nodes = buildElems(nodes, type = 'Node')
        edges = buildElems(edges, type = 'Edge')
        }

    res = shinyCyJS(c(nodes, edges), layout = layout)
    res
}


#' @import gatom
fitToBUM <- function(table, show.warnings=TRUE
) {
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }
    
    pvalsToFit <- table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
    warnWrapper(bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
    res <- bum
    res
}


lazyReadRDS <- function(name, path, envir=.GlobalEnv) {
    if (!name %in% ls(envir=envir)) {
        message(paste0("No ", name, ", loading"))
        if (startsWith(path, "http://") || startsWith(path, "https://")) {
            res <- readRDS(url(path))
        } else {
            res <- readRDS(path)
        }
        assign(name, res, envir=envir)
        message("Done")
        return(res)
    } else {
        return(get(name, envir=envir))
    }
}


normalizeName <- function(x) {
    gsub("[^a-z0-9]", "", tolower(x))
}


necessary.gene.de.fields <- c("ID", "pval", "log2FC", "baseMean")

necessary.met.de.fields <- c("ID", "pval", "log2FC")


vector2html <- function(v) {
    paste0("<ul>\n",
           paste("<li>", names(v), ": ", v, "</li>\n", collapse=""),
           "</ul>\n")
}


renderJs <- function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)

    function() {
        val <- func()
        paste0(val, ";",
               paste(sample(1:20, 10, replace=T), collapse=""))
    }
}


toJsLiteral <- function(x) {
    if (is(x, "numeric")) {
        if (x == Inf) {
            return("Infinity")
        } else {
            return("-Infinity")
        }
        return(as.character(x))
    } else if (is(x, "character")) {
        return(shQuote(x));
    } else if (is(x, "logical")) {
        return(if (x) "true" else "false")
    } else {
        stop(paste0("can't convert ", x, " to JS literal"))
    }
}


makeJsAssignments  <- function(...) {
    args <- list(...)
    values <- sapply(args, toJsLiteral)
    paste0(names(values), " = ", values, ";\n", collapse="")
}
