context("Util functions")

library(gatom)
data("org.Mm.eg.gatom.annoEx")
org.Mm.eg.gatom.anno <- org.Mm.eg.gatom.annoEx
data("gene.de.rawEx")
example.gene.de <- gene.de.rawEx
data("met.kegg.dbEx")
met.kegg.db <- met.kegg.dbEx

conf <- config::get(file=system.file('config.yml', package = 'shinyGatom'),
                    use_parent = FALSE)

test_that("lazyReadRDS works with remote files", {
    e <- new.env()
    invisible(lazyReadRDS("met.db",
                          "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds",
                          envir = e))
    expect_true(!is.null(e$met.db))
})


test_that("Node attributes are created", {
    load(system.file('testdata/moduleEx.Rda', package = 'shinyGatom'))
    vertex.table <- as_data_frame(module, what="vertices")
    res <- getJsNodeStyleAttributes(vertex.table)
    expect_true(!is.null(res$tooltip))
})


test_that("Edge attributes are created", {
    load(system.file('testdata/moduleEx.Rda', package = 'shinyGatom'))
    edge.table <- as_data_frame(module, what="edges")
    res <- getJsEdgeStyleAttributes(edge.table)
    expect_true(!is.null(res$tooltip))
})


test_that("ShinyCyJs data is created", {
    load(system.file('testdata/moduleEx.Rda', package = 'shinyGatom'))
    res <- prepareForShinyCyJS(module)
    expect_true(!is.null(res))
})


test_that("KEGG network is created with nodes as atoms", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db
    es <- makeMetabolicGraph(network=network,
                             topology=topology,
                             org.gatom.anno=org.Mm.eg.gatom.anno,
                             gene.de=example.gene.de,
                             met.db=met.db,
                             met.de=example.met.de,
                             met.to.filter=fread(system.file("mets2mask.lst",
                                                             package="gatom"))$ID)

    expect_true(!is.null(es))
})


test_that("Rhea metwork is created with nodes as metabolites", {
    network <- lazyReadRDS("rhea.network",
                           path = conf$path.to.rhea.network)
    topology <- "metabolites"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- lazyReadRDS(name = "met.rhea.db",
                          path = conf$path.to.met.rhea.db)
    es <- makeMetabolicGraph(network=network,
                             topology=topology,
                             org.gatom.anno=org.Mm.eg.gatom.anno,
                             gene.de=example.gene.de,
                             met.db=met.db,
                             met.de=example.met.de,
                             met.to.filter=fread(system.file("mets2mask.lst",
                                                             package="gatom"))$ID)

    expect_true(!is.null(es))
})


test_that("Metabolites are fitted to BUM", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst",
                                                            package="gatom"))$ID)
    k.met <- 25
    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    met.bum <- fitToBUM(table=vertex.table, k=k.met)
    expect_true(!is.null(met.bum$lambda) && !is.null(met.bum$a))
})


test_that("Genes are fitted to BUM", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst",
                                                            package="gatom"))$ID)
    k.gene <- 25
    edge.table <- data.table(as_data_frame(g, what="edges"))
    gene.bum <- fitToBUM(table=edge.table, k=k.gene)
    expect_true(!is.null(gene.bum$lambda) && !is.null(gene.bum$a))
})



test_that("Graph is scored when k.gene and k.met are not NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db

    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst",
                                                            package="gatom"))$ID)


    k.met <- 25
    k.gene <- 25
    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    met.bum <- fitToBUM(table=vertex.table, k=k.met)
    edge.table <- data.table(as_data_frame(g, what="edges"))
    gene.bum <- fitToBUM(table=edge.table, k=k.gene)
    gene.scores <- scoreNetwork(table=edge.table, k=k.gene, bum=gene.bum)
    met.scores <- scoreNetwork(table=vertex.table, k=k.met, bum=met.bum)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum,
                          vertex.threshold=met.scores$threshold, 
                          edge.threshold=gene.scores$threshold)

    expect_true(!is.null(V(sg)$score))
})



test_that("Graph is scored when k.gene is NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db

    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst",
                                                            package="gatom"))$ID)


    k.met <- 25
    k.gene <- NULL
    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    met.bum <- fitToBUM(table=vertex.table, k=k.met)
    edge.table <- data.table(as_data_frame(g, what="edges"))
    gene.bum <- fitToBUM(table=edge.table, k=k.gene)
    gene.scores <- scoreNetwork(table=edge.table, k=k.gene, bum=gene.bum)
    met.scores <- scoreNetwork(table=vertex.table, k=k.met, bum=met.bum)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum,
                          vertex.threshold=met.scores$threshold, 
                          edge.threshold=gene.scores$threshold)

    expect_true(!is.null(V(sg)$score))
})


test_that("Graph is scored when k.met is NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = conf$path.to.kegg.network)
    topology <- "atoms"
    example.met.de <- force(fread(conf$example.met.de.path))
    attr(example.met.de, "name") <- basename(conf$example.met.de.path)
    met.db <- met.kegg.db

    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst",
                                                            package="gatom"))$ID)


    k.met <- NULL
    k.gene <- 25
    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    met.bum <- fitToBUM(table=vertex.table, k=k.met)
    edge.table <- data.table(as_data_frame(g, what="edges"))
    gene.bum <- fitToBUM(table=edge.table, k=k.gene)
    gene.scores <- scoreNetwork(table=edge.table, k=k.gene, bum=gene.bum)
    met.scores <- scoreNetwork(table=vertex.table, k=k.met, bum=met.bum)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum,
                          vertex.threshold=met.scores$threshold, 
                          edge.threshold=gene.scores$threshold)

    expect_true(!is.null(V(sg)$score))
})
