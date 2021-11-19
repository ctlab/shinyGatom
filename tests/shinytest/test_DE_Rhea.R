app <- ShinyDriver$new("../../")
app$snapshotInit("test_DE_Rhea")

app$setInputs(loadExampleGeneDE = TRUE)
app$setInputs(loadExampleMetDE = TRUE)
app$setInputs(network_type = "rhea")
app$snapshot()
