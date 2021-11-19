
############################################### EDITED FROM GATOM


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
        fontSize=sapply(logPval, gatom:::getDotSize) * 25,
        bgColor=if (!is.null(attrs$log2FC)) sapply(log2FC, gatom:::getDotColor) else "white",
        borderWidth=4,
        borderColor="#eee",
        labelColor="black",
        tooltip=getJsTooltip(attrs)
    ))
}


getJsEdgeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(
        source=attrs$from,
        target=attrs$to,
        label=if (!is.null(attrs$label)) label else "",
        lineStyle="solid",
        fontSize=sapply(logPval, gatom:::getDotSize) * 25,
        lineColor=if (!is.null(attrs$log2FC)) sapply(log2FC, gatom:::getDotColor) else "grey",
        tooltip=getJsTooltip(attrs)
    ))
}


#' @import gatom
scoreGraph2 <- function(g, k.gene, k.met,
                        vertex.threshold.min=0.1,
                        edge.threshold.min=0.1,
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
        pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
        if(is.null(metabolite.bum)) {
            warnWrapper(metabolite.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
        }

        if (metabolite.bum$a > 0.5) {
            V(g)$score <- 0
            warning("Vertex scores have been assigned to 0 due to an inappropriate p-value distribution")

        } else {
            vertex.threshold <- if (k.met > length(pvalsToFit)) 1 else {
                sort(pvalsToFit)[k.met]
            }

            vertex.threshold <- min(vertex.threshold,
                                    BioNet::fdrThreshold(vertex.threshold.min, metabolite.bum))

            met.fdr <- gatom:::.reversefdrThreshold(vertex.threshold, metabolite.bum)

            gatom:::.messagef("Metabolite p-value threshold: %f", vertex.threshold)
            gatom:::.messagef("Metabolite BU alpha: %f", metabolite.bum$a)
            gatom:::.messagef("FDR for metabolites: %f", met.fdr)

            V(g)$score <- with(vertex.table,
                               (metabolite.bum$a - 1) *
                                   (log(gatom:::.replaceNA(pval, 1)) - log(vertex.threshold)))
            V(g)$score <- V(g)$score * met.score.coef
        }
    }
    else {
        V(g)$score <- 0
        V(g)$signal <- ""

    }

    if (!is.null(k.gene)) {
        pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
        if(is.null(gene.bum)){
            warnWrapper(gene.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
        }

        if(gene.bum$a > 0.5) {
            E(g)$score <- 0
            warning("Edge scores have been assigned to 0 due to an inappropriate p-value distribution")

        } else {
            edge.threshold <- if (k.gene > length(pvalsToFit)) 1 else {
                sort(pvalsToFit)[k.gene]
            }

            edge.threshold <- min(edge.threshold,
                                  BioNet::fdrThreshold(edge.threshold.min, gene.bum))

            gene.fdr <- gatom:::.reversefdrThreshold(edge.threshold, gene.bum)

            gatom:::.messagef("Gene p-value threshold: %f", edge.threshold)
            gatom:::.messagef("Gene BU alpha: %f", gene.bum$a)
            gatom:::.messagef("FDR for genes: %f", gene.fdr)

            E(g)$score <- with(edge.table,
                               (gene.bum$a - 1) *
                                   (log(gatom:::.replaceNA(pval, 1)) - log(edge.threshold)))
        }
    }
    else {
        E(g)$score <- 0
        E(g)$signal <- ""
    }
    g
    if (raw) {
        return(g)
    }

    res <- normalize_sgmwcs_instance(g,
                                     nodes.weight.column = "score",
                                     edges.weight.column = "score",
                                     nodes.group.by = "signal",
                                     edges.group.by = "signal",
                                     group.only.positive = TRUE)
    res
}


############################################### NEW FUNCTIONS

prepareForshinyCyJS <- function(module){
    if (is.null(module)) {
        return(NULL)
    }
    
    vertex.table <- as_data_frame(module, what="vertices")
    edge.table <- as_data_frame(module, what="edges")
    
    nodes <- getJsNodeStyleAttributes(vertex.table)
    edges <- getJsEdgeStyleAttributes(edge.table)
    
    nodes = buildElems(nodes, type = 'Node')
    edges = buildElems(edges, type = 'Edge')
    
    res = shinyCyJS(c(nodes, edges))
    res
}


#' @import gatom
fitGenesToBUM <- function(g, 
                          k.gene,
                          show.warnings=TRUE
) {
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }

    edge.table <- data.table(as_data_frame(g, what="edges"))

    if (!is.null(k.gene)) {
        pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
        warnWrapper(gene.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
    }
    res <- gene.bum
    res
}


#' @import gatom
fitMetsToBUM <- function(g, 
                         k.met,
                         show.warnings=TRUE
) {
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }

    vertex.table <- data.table(as_data_frame(g, what="vertices"))

    if (!is.null(k.met)) {
        pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]
        warnWrapper(metabolite.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
    }
    res <- metabolite.bum
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


############################################### OLD FUNCTIONS


normalizeName <- function(x) {
    gsub("[^a-z0-9]", "", tolower(x))
}


read.table.smart <- function(path, ...) {
    fields <- list(...)
    conn <- file(path)
    header <- readLines(conn, n=1)
    close(conn)

    seps <- c("\t", " ", ",", ";")
    sep <- seps[which.max(table(unlist(strsplit(header, "")))[seps])]

    res <- read.table(path, sep=sep, header=T, stringsAsFactors=F, check.names=F, quote='"')
    res <- as.data.table(res, keep.rownames=is.character(attr(res, "row.names")))

    oldnames <- character(0)
    newnames <- character(0)

    for (field in names(fields)) {
        if (field %in% colnames(res)) {
            next
        }

        z <- na.omit(
            match(
                normalizeName(c(field, fields[[field]])),
                normalizeName(colnames(res))))
        if (length(z) == 0) {
            next
        }

        oldnames <- c(oldnames, colnames(res)[z[1]])
        newnames <- c(newnames, field)
    }

    logdebug("smart renaming")
    logdebug("from: %s", paste0(oldnames, collapse=" "))
    logdebug("to: %s", paste0(newnames, collapse=" "))
    setnames(res, oldnames, newnames)
    res
}


read.table.smart.de <- function(path, ID=ID) {
    read.table.smart(path, ID=ID, pval=c("pvalue"), log2FC=c("log2foldchange", "logfc"), name.orig=c("name"))
}


read.table.smart.de.gene <- function(path, idsList) {
    res <- read.table.smart.de(path, ID="ID")
    idColumn <- findIdColumn(res, idsList)
    if (idColumn$matchRatio < 0.1) {
        z <- na.omit(
            match(
                normalizeName(c("ID", "gene id", "gene", "entrez", "", "rn", "symbol")),
                normalizeName(colnames(res))))
        if (length(z) == 0) {
            setnames(res, colnames(res)[1], "ID")
        } else {
            setnames(res, colnames(res)[z[1]], "ID")
        }

    } else {
        if (idColumn$column != "ID" && "ID" %in% colnames(res)) {
            setnames(res, "ID", "ID.old")
        }
        setnames(res, idColumn$column, "ID")
    }
    res
}

read.table.smart.de.met <- function(path) {
    res <- read.table.smart.de(path, ID=c("metabolite", "kegg", "hmdb", "", "rn"))
    if (!"ID" %in% colnames(res)) {
        setnames(res, colnames(res)[1], "ID")
    }
    res
}


necessary.de.fields <- c("ID", "pval")


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


.intersectionSize <- function(...) { length(intersect(...))}


findIdColumn <- function(de, idsList,
                         sample.size=1000,
                         match.threshold=0.6,
                         remove.ensembl.revisions=TRUE) {
    # first looking for column with base IDs
    de.sample <- if (nrow(de) < sample.size) {
        copy(de)
    } else {
        de[sample(seq_len(nrow(de)), sample.size), ]
    }
    columnSamples <- lapply(de.sample, as.character)


    if (remove.ensembl.revisions) {
        columnSamples <- lapply(columnSamples, gsub,
                                pattern="(ENS\\w*\\d*)\\.\\d*",
                                replacement="\\1")
    }

    ss <- sapply(columnSamples,
                 .intersectionSize, idsList[[1]])

    if (max(ss) / nrow(de.sample) >= match.threshold) {
        # we found a good column with base IDs
        return(list(column=colnames(de)[which.max(ss)],
                    type=names(idsList)[1],
                    matchRatio=max(ss) / nrow(de.sample)))
    }

    z <- .pairwiseCompare(.intersectionSize,
                          columnSamples,
                          idsList)

    bestMatch <- which(z == max(z), arr.ind = TRUE)[1,]
    return(list(column=colnames(de)[bestMatch["row"]],
                type=names(idsList)[bestMatch["col"]],
                matchRatio=max(z) / nrow(de.sample)))
}


.pairwiseCompare <- function (FUN, list1, list2 = list1, ...)
{
    additionalArguments <- list(...)
    f1 <- function(...) {
        mapply(FUN = function(x, y) {
            do.call(FUN, c(list(list1[[x]], list2[[y]]), additionalArguments))
        }, ...)
    }
    z <- outer(seq_along(list1), seq_along(list2), FUN = f1)
    rownames(z) <- names(list1)
    colnames(z) <- names(list2)
    z
}
