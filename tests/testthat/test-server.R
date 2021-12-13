library(shiny)

# devtools::load_all()
# debug(testServer)
# undebug(testServer)


## Data loads & annotates with KEGG network
context("Data loads & annotates with KEGG network")
testServer(ShinyGATOM(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "kegg")

    expect_true(metIdsType() == "HMDB")
})


## Data loads & annotates with Rhea network
context("Data loads & annotates with Rhea network")
testServer(ShinyGATOM(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "rhea")

    expect_true(metIdsType() == "HMDB")
})


## Network is created with KEGG network and nodes as atoms
context("Network is created with KEGG network and nodes as atoms")
testServer(ShinyGATOM(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "kegg")
    session$setInputs(organism = "mmu")
    session$setInputs(nodesAs = "atoms")
    session$setInputs(autoFindModule = TRUE)
    session$setInputs(preprocess = TRUE)
    session$setInputs(runStep1 = "click")
    g <- gInput()
    expect_true(!is.null(V(g)))
})


## Network is created with Rhea network and nodes as metabolites
context("Network is created with Rhea network and nodes as metabolites")
testServer(ShinyGATOM(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "rhea")
    session$setInputs(organism = "mmu")
    session$setInputs(nodesAs = "metabolites")
    session$setInputs(autoFindModule = TRUE)
    session$setInputs(preprocess = TRUE)
    session$setInputs(runStep1 = "click")
    g <- gInput()
    expect_true(!is.null(V(g)))
})


# ## Module is found
# context("Module is found")
# testServer(ShinyGATOM(), {
#     session$setInputs(loadExampleGeneDE = TRUE)
#     session$setInputs(loadExampleMetDE = TRUE)
#     session$setInputs(network = "kegg")
#     session$setInputs(organism = "mmu")
#     session$setInputs(nodesAs = "atoms")
#     session$setInputs(autoFindModule = TRUE)
#     session$setInputs(preprocess = TRUE)
#     
#     session$setInputs(kgene = 25)
#     session$setInputs(nullkgene = FALSE)
#     session$setInputs(kmet = 25)
#     session$setInputs(nullkmet = FALSE)
#     
#     session$setInputs(runStep1 = "click")
#     
#     session$setInputs(addHighlyExpressedEdges = FALSE)
#     session$setInputs(solveToOptimality = FALSE)
#     session$setInputs(metaboliteActions = "NoAction")
#     session$setInputs(find = TRUE)
#     
#     session$setInputs(runAll = "click")
#     
#     module <- moduleInput()
#     expect_true(!is.null(V(module)))
# })


# context("")