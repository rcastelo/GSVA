test_sparseMethods <- function(){
    message("Running unit tests for sparse methods.")
    
    set.seed(123)
    
    m <- matrix(sample.int(10, 25, T), 10, 10)
    colnames(m) <- paste0("cell_", 1:10)
    rownames(m) <- paste0("gene_", 1:10)
    
    gene.sets <- list("my_list1"= paste0("gene_", 1:2),
                      "my_list2"= paste0("gene_", 3:4))
    
    M <- as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix")

    mg <- gsva(gsvaParam(m, gene.sets), verbose=FALSE)
    Mg <- gsva(gsvaParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(mg, Mg)
    
    mp <- gsva(plageParam(m, gene.sets), verbose=FALSE)
    Mp <- gsva(plageParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(mp, Mp)
    
    mz <- gsva(zscoreParam(m, gene.sets), verbose=FALSE)
    Mz <- gsva(zscoreParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(mz, Mz)
    
    ms <- gsva(ssgseaParam(m, gene.sets), verbose=FALSE)
    Ms <- gsva(ssgseaParam(M, gene.sets), verbose=FALSE)
    checkEqualsNumeric(ms, Ms)
}
