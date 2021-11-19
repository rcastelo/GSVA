test_delayed <- function(){
  
  # set.seed(1001)
  m <- matrix(runif(10000), 1000, 100)
  rownames(m) <- paste0("gene_", 1:1000)
  colnames(m) <- paste0("cell_", 1:100)
  
  gene.sets <- list("my_list1"= paste0("gene_", 1:100),
                    "my_list2"= paste0("gene_", 401:500))
  
  h5 <- as(m, "HDF5Array")
  
  x1 <- gsva(m, gene.sets, method="plage")
  x2 <- gsva(h5, gene.sets, method="plage")
  
  ## SVD sign is arbitrary and implementation-dependent
  checkEqualsNumeric(abs(as.vector(t(x1))), abs(as.vector(t(x2))))
  
  x1 <- gsva(m, gene.sets, method="zscore")
  x2 <- gsva(h5, gene.sets, method="zscore")
  
  checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
  
  x1 <- gsva(m, gene.sets, method="ssgsea")
  x2 <- gsva(h5, gene.sets, method="ssgsea")
  
  checkEqualsNumeric(as.vector(t(x1)), as.vector(t(x2)))
}
