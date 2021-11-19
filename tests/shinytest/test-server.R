testServer(expr = {
    # Set the `size` slider and check the output
    session$loadExampleGeneDE(value=TRUE)
    session$loadExampleMetDE(value=TRUE)
    session$network_type(selected="kegg")
    
    
    session$setInputs(size = 12)
    expect_equal(output$sequence, "1 10 11 12 2 3 4 5 6 7 8 9")
    
})