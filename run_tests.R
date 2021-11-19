library(testthat)
library(shinytest)

test_that("KEGG network with nodes as metabolites works", {
    # Use compareImages=FALSE because the expected image screenshots were created
    # on a Mac, and they will differ from screenshots taken on the CI platform,
    # which runs on Linux.
    expect_pass(testApp(".", testnames="KEGG_module_testing", compareImages = FALSE))
})

test_that("KEGG network with nodes as atoms works", {
    # Use compareImages=FALSE because the expected image screenshots were created
    # on a Mac, and they will differ from screenshots taken on the CI platform,
    # which runs on Linux.
    expect_pass(testApp(".", testnames="test_DE_KEGG", compareImages = FALSE))
})


test_that("KEGG module works", {
    # Use compareImages=FALSE because the expected image screenshots were created
    # on a Mac, and they will differ from screenshots taken on the CI platform,
    # which runs on Linux.
    expect_pass(testApp(".", testnames="KEGG_module_testing", compareImages = FALSE))
})
