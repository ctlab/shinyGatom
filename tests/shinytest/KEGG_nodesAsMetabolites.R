app <- ShinyDriver$new("../../")
app$snapshotInit("KEGG_nodesAsMetabolites")

app$setInputs(loadExampleGeneDE = TRUE)
app$setInputs(loadExampleMetDE = TRUE)
app$setInputs(nodesAs = "metabolites")
app$snapshot()
app$snapshot()
