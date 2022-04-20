library(shiny)

# devtools::load_all()
# debug(testServer)
# undebug(testServer)
# devtools::test_active_file()

# Sys.setenv("R_CONFIG_ACTIVE"="local")
# Sys.setenv("GATOM_LOCAL_DATA"="/home/masha/Research/RCode/RShiny/data")
# Sys.setenv("CPLEX_HOME"="/home/masha/cplex")

conf <- config::get(file=system.file('config.yml', package = 'shinyGatom'),
                    use_parent = FALSE)


## Data is loaded & annotated with KEGG network
context("Data loads & annotates with KEGG network")
testServer(shinyGatom(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "kegg")
    
    expect_true(metIdsType() == "HMDB")
})


## Data is read & annotated with KEGG network
context("Data is read & annotated with KEGG network")
testServer(shinyGatom(), {
    session$setInputs(loadExampleGeneDE = FALSE)
    session$setInputs(loadExampleMetDE = FALSE)
    session$setInputs(loadExampleLipidDE = FALSE)
    
    session$setInputs(network = "kegg")
    session$setInputs(organism = "mmu")
    session$setInputs(nodesAs = "atoms")
    
    session$setInputs(geneDE=data.frame(datapath=conf$example.gene.de.path,
                                        name="Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    session$setInputs(metDE=data.frame(datapath=conf$example.met.de.path,
                                       name="Ctrl.vs.MandLPSandIFNg.met.de.tsv"))
    
    expect_equal(metIdsType(), "HMDB")
})


## Data is loaded & annotated with Rhea network
context("Data is loaded & annotated with Rhea network")
testServer(shinyGatom(), {
    session$setInputs(loadExampleGeneDE = TRUE)
    session$setInputs(loadExampleMetDE = TRUE)
    session$setInputs(loadExampleLipidDE = FALSE)
    session$setInputs(network = "rhea")
    
    expect_true(metIdsType() == "HMDB")
})


## Network is created with KEGG network and nodes as atoms
context("Network is created with KEGG network and nodes as atoms")
testServer(shinyGatom(), {
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
testServer(shinyGatom(), {
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


context("Network is created with KEGG network and nodes as metabolites")
testServer(shinyGatom(), {
    session$setInputs(loadExampleGeneDE = FALSE)
    session$setInputs(loadExampleMetDE = FALSE)
    session$setInputs(loadExampleLipidDE = FALSE)
    
    session$setInputs(network = "kegg")
    session$setInputs(organism = "hsa")
    session$setInputs(nodesAs = "metabolites")
    
    session$setInputs(geneDE=data.frame(datapath="https://artyomovlab.wustl.edu/publications/supp_materials/GAM/MCF10A.Ctrl.vs.2DG.gene.de.tsv",
                                        name="MCF10A.Ctrl.vs.2DG.gene.de.tsv"))
    
    session$setInputs(preprocess = TRUE)
    session$setInputs(runStep1 = "click")
    
    g <- gInput()
    expect_true(!is.null(V(g)))
})



# ## Module is found
# context("Module is found")
# testServer(shinyGatom(), {
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
