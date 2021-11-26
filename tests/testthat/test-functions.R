context("Util functions")

# devtools::test()
devtools::load_all()

test_that("lazyReadRDS works with remote files", {
    e <- new.env()
    invisible(lazyReadRDS("met.db", "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds", envir = e))
    expect_true(!is.null(e$met.db))
})


test_that("Node attributes are created", {
    load("/home/masha/Research/RCode/RShiny/example_module2.Rda")
    vertex.table <- as_data_frame(module, what="vertices")
    res <- getJsNodeStyleAttributes(vertex.table)
    expect_true(!is.null(res$tooltip))
})


test_that("Edge attributes are created", {
    load("/home/masha/Research/RCode/RShiny/example_module2.Rda")
    edge.table <- as_data_frame(module, what="edges")
    res <- getJsEdgeStyleAttributes(edge.table)
    expect_true(!is.null(res$tooltip))
})


test_that("ShinyCyJs data is created", {
    load("/home/masha/Research/RCode/RShiny/example_module2.Rda")
    res <- prepareForShinyCyJS(module)
    expect_true(!is.null(res))
})


test_that("KEGG network is created with nodes as atoms", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    es <- makeMetabolicGraph(network=network,
                             topology=topology,
                             org.gatom.anno=org.Mm.eg.gatom.anno,
                             gene.de=example.gene.de,
                             met.db=met.db,
                             met.de=example.met.de,
                             met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    expect_true(!is.null(es))
})


test_that("Rhea metwork is created with nodes as metabolites", {
    network <- lazyReadRDS("rhea.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.rhea.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "metabolites"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.rhea.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.rhea.db.rds")
    es <- makeMetabolicGraph(network=network,
                             topology=topology,
                             org.gatom.anno=org.Mm.eg.gatom.anno,
                             gene.de=example.gene.de,
                             met.db=met.db,
                             met.de=example.met.de,
                             met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    expect_true(!is.null(es))
})


test_that("Metabolites are fitted to BUM", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    met.bum <- fitMetsToBUM(g=g, k.met=25)
    expect_true(!is.null(met.bum$lambda) && !is.null(met.bum$a))
})


test_that("Genes are fitted to BUM", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    gene.bum <- fitGenesToBUM(g=g, k.gene=25)
    expect_true(!is.null(gene.bum$lambda) && !is.null(gene.bum$a))
})



test_that("Graph is scored when k.gene and k.met are not NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    
    k.met <- 25
    k.gene <- 25
    met.bum <- fitMetsToBUM(g=g, k.met=k.met)
    gene.bum <- fitGenesToBUM(g=g, k.gene=k.gene)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum)
    
    expect_true(!is.null(V(sg)$score))
})



test_that("Graph is scored when k.gene is NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    
    k.met <- 25
    k.gene <- NULL
    met.bum <- fitMetsToBUM(g=g, k.met=k.met)
    gene.bum <- fitGenesToBUM(g=g, k.gene=k.gene)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum)
    
    expect_true(!is.null(V(sg)$score))
})


test_that("Graph is scored when k.met is NULL", {
    network <- lazyReadRDS("kegg.network",
                           path = "/home/masha/Research/RCode/RShiny/data/network.kegg.rds")
    org.Mm.eg.gatom.anno <- lazyReadRDS("org.Mm.eg.gatom.anno",
                                        "/home/masha/Research/RCode/RShiny/data/org.Mm.eg.gatom.anno.rds")
    topology <- "atoms"
    example.gene.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"))
    attr(example.gene.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv")
    example.met.de <- force(fread("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    attr(example.met.de, "name") <- basename("/home/masha/Research/RCode/RShiny/data/Ctrl.vs.MandLPSandIFNg.met.de.tsv")
    met.db <- lazyReadRDS(name = "met.kegg.db",
                          path = "/home/masha/Research/RCode/RShiny/data/met.kegg.db.rds")
    
    g <- makeMetabolicGraph(network=network,
                            topology=topology,
                            org.gatom.anno=org.Mm.eg.gatom.anno,
                            gene.de=example.gene.de,
                            met.db=met.db,
                            met.de=example.met.de,
                            met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID)
    
    
    k.met <- NULL
    k.gene <- 25
    met.bum <- fitMetsToBUM(g=g, k.met=k.met)
    gene.bum <- fitGenesToBUM(g=g, k.gene=k.gene)
    
    sg <- scoreGraphShiny(g=g, k.gene=k.gene, k.met=k.met,
                          metabolite.bum=met.bum, gene.bum=gene.bum)
    
    expect_true(!is.null(V(sg)$score))
})
