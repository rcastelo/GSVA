source(system.file("unitTests", "test_helpers.R", package = "GSVA"))

test_walkBins <- function(){
    message("Running unit tests for the GPU walk size ladder")

    p <- 4300L                      ## genes; must exceed the top capacity 4096
    n <- 12L                        ## samples
    edges <- c(31L, 32L, 33L,       ## category 0: capacity 32
               63L, 64L, 65L,       ## category 1: capacity 64
               127L, 128L, 129L,    ## category 2: capacity 128
               159L, 160L, 161L,    ## category 3: capacity 160
               191L, 192L, 193L,    ## category 4: capacity 192
               255L, 256L, 257L,    ## category 5: capacity 256
               383L, 384L, 385L,    ## category 6: capacity 384
               767L, 768L, 769L,    ## category 7: capacity 768
               1535L, 1536L, 1537L, ## category 8: capacity 1536
               2047L, 2048L, 2049L, ## category 9, lower half
               4095L, 4096L)        ## category 10: capacity 4096

    set.seed(42)
    x <- matrix(rnorm(p * n), nrow = p, ncol = n,
                dimnames = list(paste0("g", seq_len(p)), paste0("s", seq_len(n))))
    gsets <- lapply(edges, function(k) paste0("g", seq_len(k)))
    names(gsets) <- paste0("gs", edges)

    gsvaPar <- gsvaParam(x, gsets, verbose = FALSE)
    cpu.es <- gsva(gsvaPar, verbose = FALSE, device = "cpu")
    .run_with_gpu({
        gpu.es <- gsva(gsvaPar, verbose = FALSE, device = "gpu")
        .check_gpu_equiv(cpu.es, gpu.es)
    })
}

test_walkBinsSparse <- function(){
    message("Running unit tests for the GPU sparse walk size ladder")

    suppressPackageStartupMessages(library(Matrix))
    p <- 4300L                      ## genes; must exceed the top capacity 4096
    n <- 12L                        ## samples
    edges <- c(31L, 32L, 33L,       ## capacity 32
               63L, 64L, 65L,       ## capacity 64
               127L, 128L, 129L,    ## capacity 128
               159L, 160L, 161L,    ## capacity 160
               191L, 192L, 193L,    ## capacity 192
               255L, 256L, 257L,    ## capacity 256
               383L, 384L, 385L,    ## capacity 384
               767L, 768L, 769L,    ## capacity 768
               1535L, 1536L, 1537L, ## capacity 1536
               2047L, 2048L, 2049L, ## capacity 2048
               4095L, 4096L)        ## capacity 4096

    set.seed(42)
    x <- matrix(ifelse(runif(p * n) < 0.2,
                       round(runif(p * n, 1, 500)), 0),
                nrow = p, ncol = n,
                dimnames = list(paste0("g", seq_len(p)), paste0("s", seq_len(n))))
    x <- Matrix(x, sparse = TRUE)
    gsets <- lapply(edges, function(k) paste0("g", seq_len(k)))
    names(gsets) <- paste0("gs", edges)

    gsvaPar <- gsvaParam(x, gsets, verbose = FALSE, sparse = TRUE)
    cpu.es <- gsva(gsvaPar, verbose = FALSE, device = "cpu")
    .run_with_gpu({
        gpu.es <- gsva(gsvaPar, verbose = FALSE, device = "gpu")
        .check_gpu_equiv(cpu.es, gpu.es)
    })
}
