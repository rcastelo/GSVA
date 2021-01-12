test_delayedMethods <- function(){
  
  m <- matrix(runif(100000), 10000, 1000)
  rownames(m) <- paste0("gene_", 1:10000)
  colnames(m) <- paste0("cell_", 1:1000)
  
  gene.sets <- list("my_list1"= paste0("gene_", 1:1000),
                    "my_list2"= paste0("gene_", 4001:5000))
  
  h5 <- as(m, "HDF5Array")
  
  x1 <- gsva(m, gene.sets, method="plage")
  x2 <- gsva(h5, gene.sets, method="plage")
  
  checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
  
  x1 <- gsva(m, gene.sets, method="zscore")
  x2 <- gsva(h5, gene.sets, method="zscore")
  
  checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
  
  x1 <- gsva(m, gene.sets, method="ssgsea")
  x2 <- gsva(h5, gene.sets, method="ssgsea")
  
  checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
}




