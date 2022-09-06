test_sparseMethods <- function(){
  
  set.seed(123)
  
  m <- matrix(sample.int(10, 25, T), 10, 10)
  colnames(m) <- paste0("cell_", 1:10)
  rownames(m) <- paste0("gene_", 1:10)
  
  gene.sets <- list("my_list1"= paste0("gene_", 1:2),
                    "my_list2"= paste0("gene_", 3:4))
  
  M <- as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  
  checkEqualsNumeric(gsva(m, gene.sets), gsva(M, gene.sets))
  checkEqualsNumeric(gsva(m, gene.sets, method="plage"), gsva(M, gene.sets, method="plage"))
  checkEqualsNumeric(gsva(m, gene.sets, method="zscore"), gsva(M, gene.sets, method="zscore"))
  checkEqualsNumeric(gsva(m, gene.sets, method="ssgsea"), gsva(M, gene.sets, method="ssgsea"))
}
