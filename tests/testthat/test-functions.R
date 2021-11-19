context("Util functions")

test_that("lazyeReadRDS works with remote files", {
    e <- new.env()
    invisible(lazyReadRDS("met.db", "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds", envir = e))
    expect_true(!is.null(e$met.db))
})
