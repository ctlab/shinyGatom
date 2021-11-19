app <- ShinyDriver$new("../../")
app$snapshotInit("test_DE_KEGG")

app$setInputs(loadExampleGeneDE = TRUE)
app$setInputs(loadExampleMetDE = TRUE)
app$snapshot()
