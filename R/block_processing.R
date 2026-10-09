## This file contains the internal functions that are used to process the rows
## or the columns of a matrix in blocks, stored in main memory or in an on-disk
## data structure such as an HDF5Matrix object, opening parallelism through a
## BiocParallelParam object when different from NULL, and reporting progress
## using the 'cli' package when possible. The block size is automatically
## determined based on the number of workers and the maximum available memory,
## if specified. The functions in this file are not exported and are not
## intended to be used directly by the users of the package. They are only used
## internally by the functions that need to process the rows or the columns of
## a matrix in blocks.

## fraction of the maximum main memory available to the allocations of R, the
## rest being taken by the memory that R releases but that the process keeps,
## and by allocations outside R, e.g., by the HDF5 library; in the steps of
## GSVA on single-cell data, the peak of memory allocated by R was 70% to 80%
## of the peak of memory resident in the process
.mem_fraction_R <- 0.7

## memory taken by the steps of GSVA to process each row or column of the
## input data 'X', relative to its size in 'X' ('workfactor'), and by the
## output of each row or column, relative to that size ('outfactor') plus a
## number of bytes ('outextra'), with 'ngs' gene sets in the steps giving
## enrichment scores, which are dense. the
## factors were measured on single-cell data and include a margin; the output
## of normalizing rows is double, and of ranking columns integer. instead of
## 'X', whether its values are sparse or integer can be given directly. the
## optional factor 'assembly' gives the memory taken by the output in main
## memory, relative to its size, while it is produced, see .fixed_mem(). the
## row normalization with CLR ('clr=TRUE') reads its input once through blocks
## of columns, see .rowStats(), and gives its output in one piece, which takes
## twice its size when the input is a dgCMatrix object ('dgc=TRUE') converted
## into an SVT_SparseMatrix object, see .rownorm_clr(); when the input is
## stored on disk, its output is a delayed operation that takes no memory
#' @importFrom BiocGenerics type
#' @importFrom S4Arrays is_sparse
.step_mem_factors <- function(step=c("rownorm", "rowstats", "colranks",
                                     "scores", "average", "plage", "zscore",
                                     "ssgsea", "ssgseascores"), X=NULL, ngs=0,
                              sparse=is_sparse(X),
                              int=(type(X) == "integer"), clr=FALSE,
                              dgc=is(X, "dgCMatrix")) {
    step <- match.arg(step)
    if (step == "rownorm" && clr)
        return(list(workfactor=if (sparse) 2 else 5,
                    outfactor=if (!int) 1 else if (sparse) 1.5 else 2,
                    outextra=0, assembly=if (dgc) 2 else 1))
    switch(step,
           rownorm=list(workfactor=if (sparse) 6 else 4,
                        outfactor=if (!int) 1 else if (sparse) 1.5 else 2,
                        outextra=0),
           ## the row statistics of each block of columns, four or five
           ## doubles per row, are small and reduced at the end, see
           ## .rowStats(); dense blocks take a copy of the block, the mask of
           ## missing values, the logarithm of their values and temporary
           ## matrices, about five times their size, as measured
           rowstats=list(workfactor=if (sparse) 2 else 5, outfactor=0,
                         outextra=0),
           ## the ranks of a dgCMatrix object are double, as large as it,
           ## and binding their blocks of columns takes twice their size,
           ## see .cbind_blocks()
           colranks=list(workfactor=if (sparse) 7 else 2,
                         outfactor=if (int || (sparse && dgc)) 1
                                   else if (sparse) 0.7 else 0.5,
                         outextra=0, assembly=if (sparse && dgc) 2 else NULL),
           scores=list(workfactor=2, outfactor=0, outextra=8 * ngs),
           ## the average method gives its scores directly, while PLAGE,
           ## z-score and ssGSEA give them from an intermediate matrix as large
           ## as their input, of scaled rows or column ranks
           average=list(workfactor=2, outfactor=0, outextra=8 * ngs),
           plage=, zscore=,
           ssgsea=list(workfactor=if (step == "ssgsea") 3 else 2, outfactor=1,
                       outextra=8 * ngs),
           ## the scores of ssGSEA, from its ranks, which take from four to
           ## six times the size of a block of integer ranks: the block read
           ## from disk, its integer copy, its ranks to the power of alpha and
           ## the temporary vectors of its columns, as measured on single-cell
           ## data, the larger factors with smaller blocks
           ssgseascores=list(workfactor=6, outfactor=0, outextra=8 * ngs))
}

## memory that remains allocated while a matrix with 'nunits' rows (whdim=1)
## or columns (whdim=2) of 'unitbytes' bytes each, sparse or not, is processed
## in blocks: the matrix itself, when it is in main memory, and the output,
## see .step_mem_factors(). when the matrix is stored on disk, the output
## proportional to it is written to disk block by block. otherwise, the output
## is assembled from the blocks at the end, which takes twice its size, except
## when binding blocks of columns of sparse data, which reuses the memory of
## the blocks, or as given by 'assembly'. the output proportional to the
## matrix takes 'outunitbytes' per row or column, when its size differs from
## the one of the input, such as dense output from sparse input. the output
## given by 'outextra' is always assembled in memory
.fixed_mem <- function(nunits, unitbytes, whdim, sparse, inmemory, outfactor,
                       outextra, assembly=NULL, outunitbytes=unitbytes) {
    fixed <- 2 * outextra * nunits
    if (inmemory) {
        if (is.null(assembly))
            assembly <- if (whdim == 2L && sparse) 1 else 2
        fixed <- fixed + nunits * unitbytes +
                 assembly * outfactor * outunitbytes * nunits
    }

    fixed
}

## number of rows (whdim=1) or columns (whdim=2) of the blocks in which the
## matrix 'X' can be processed by 'nworkers' workers within the maximum main
## memory 'maxmem', where 'workfactor', 'outfactor' and 'outextra' are given
## by .step_mem_factors(). the memory available to the blocks of the workers
## is what is left by the memory that remains allocated, see .fixed_mem(),
## from the fraction .mem_fraction_R, available to the allocations of R, of
## the memory that the R processes of more than one worker, of 'workermem'
## bytes each, see .worker_mem(), leave from 'maxmem', when they leave enough
## for the memory that remains allocated; 'workermem=0' when the workers do
## not share 'maxmem', such as the jobs of gsvaMap(). the size of a row or column of 'X' stored on disk is the one of its
## dense form, which overestimates the size of sparse data. with 'dense=TRUE',
## the calculations convert each block of a sparse 'X' in main memory into a
## dense matrix, whose rows or columns take their working memory and output
## by the size of their dense form, while 'X' itself takes its sparse size.
## 'heldmem' bytes remain allocated by other objects, such as the input of
## previous steps still in main memory, see .held_mem(). the blocks of all
## workers together are not smaller than the automatic block size of the
## DelayedArray package, see getAutoBlockSize(), because the overhead of
## processing many smaller blocks makes calculations too slow, so that a
## smaller 'maxmem' is not honored, and
## they have at most as many rows or columns as needed to give one block to
## each worker, and as the DelayedArray package supports
#' @importFrom BiocGenerics type
#' @importFrom utils object.size
#' @importFrom DelayedArray getAutoBlockLength
#' @importFrom S4Arrays is_sparse
.units_per_block <- function(X, whdim, nworkers, maxmem, workfactor=2,
                             outfactor=1, outextra=0,
                             workermem=.worker_mem(), dense=FALSE,
                             heldmem=.held_mem()) {
    nunits <- dim(X)[whdim]
    inmemory <- !is(X, "DelayedArray")
    eltbytes <- if (type(X) == "integer") 4 else 8
    denseunitbytes <- as.numeric(dim(X)[-whdim]) * eltbytes
    if (inmemory)
        inunitbytes <- as.numeric(object.size(X)) / max(1, nunits)
    else
        inunitbytes <- denseunitbytes
    sparse <- is_sparse(X) && !dense
    unitbytes <- if (dense) denseunitbytes else inunitbytes
    ## binding blocks of columns of a dgCMatrix object takes twice the size
    ## of the result, as measured, unlike binding SVT_SparseMatrix objects
    assembly <- if (whdim == 2L && sparse && is(X, "dgCMatrix")) 2 else NULL
    fixed <- .fixed_mem(nunits, inunitbytes, whdim, sparse, inmemory,
                        outfactor, outextra, assembly,
                        outunitbytes=unitbytes) + heldmem
    ## with more than one worker, their R processes take memory from 'maxmem',
    ## unless they leave no memory for the blocks, when 'maxmem' cannot be
    ## honored anyway, and smaller blocks would only make calculations slower
    if (nworkers > 1 &&
        .mem_fraction_R * (maxmem - nworkers * workermem) > fixed)
        maxmem <- maxmem - nworkers * workermem
    avail <- (.mem_fraction_R * maxmem - fixed) / nworkers
    nperblock <- floor(avail / ((workfactor + outfactor) * unitbytes +
                                outextra))
    otherdim <- max(1, as.numeric(dim(X)[-whdim]))
    nperblock <- max(nperblock, floor(getAutoBlockLength(type(X)) / nworkers /
                                      otherdim))
    ## the DelayedArray package does not support blocks with more than
    ## .Machine$integer.max values
    nperblock <- min(ceiling(nunits / nworkers),
                     floor(.Machine$integer.max / otherdim), nperblock)

    as.integer(max(1, nperblock))
}

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R; with a
## finite maximum main memory 'maxmem', blocks have as many rows as fit in it,
## see .units_per_block(), and otherwise the default block length of the
## DelayedArray package, when there are several workers or 'X' is stored on
## disk, or form a single block
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.rowgridsize <- function(X, nworkers=1, maxmem=Inf, workfactor=2, outfactor=1,
                         outextra=0, workermem=.worker_mem(), dense=FALSE,
                         heldmem=.held_mem()) {
  grid <- DummyArrayGrid(dim(X))
  if (is.finite(maxmem)) {
      nrowblock <- .units_per_block(X, 1L, nworkers, maxmem, workfactor,
                                    outfactor, outextra, workermem,
                                    dense=dense, heldmem=heldmem)
      grid <- rowAutoGrid(X, nrow=.align_to_chunks(nrowblock, X, 1L))
  } else if (nworkers > 1 || is(X, "DelayedMatrix")) {
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers and, in any case, it cannot exceed
      ## .Machine$integer.max
      max.block.length <- min(.Machine$integer.max,
                              getAutoBlockLength(type(X)) / nworkers)
      expected.block.length <- max(1, ceiling(nrow(X) / nworkers) * as.numeric(ncol(X)))
      block.length <- min(max.block.length, expected.block.length)
      ## number of rows per block, as calculated by rowAutoGrid()
      nrowblock <- min(max(1, floor(block.length / max(1, ncol(X)))), nrow(X))
      grid <- rowAutoGrid(X, nrow=.align_to_chunks(nrowblock, X, 1L))
  }
  grid
}

## adapted from .define_multiworker_grid() in beachmat/R/colBlockApply.R; with a
## finite maximum main memory 'maxmem', blocks have as many columns as fit in
## it, see .units_per_block(), and otherwise the default block length of the
## DelayedArray package, when there are several workers or 'X' is stored on
## disk, or form a single block
#' @importFrom BiocGenerics type
#' @importFrom DelayedArray rowAutoGrid colAutoGrid getAutoBlockLength
.colgridsize <- function(X, nworkers=1, maxmem=Inf, workfactor=2, outfactor=1,
                         outextra=0, workermem=.worker_mem(), dense=FALSE,
                         heldmem=.held_mem()) {
  grid <- DummyArrayGrid(dim(X))
  if (is.finite(maxmem)) {
      ncolblock <- .units_per_block(X, 2L, nworkers, maxmem, workfactor,
                                    outfactor, outextra, workermem,
                                    dense=dense, heldmem=heldmem)
      grid <- colAutoGrid(X, ncol=.align_to_chunks(ncolblock, X, 2L))
  } else if (nworkers > 1 || is(X, "DelayedMatrix")) {
      ## assuming all workers share memory, the maximum block length has to reduce
      ## by the number of workers and, in any case, it cannot exceed
      ## .Machine$integer.max
      max.block.length <- min(.Machine$integer.max,
                              getAutoBlockLength(type(X)) / nworkers)
      expected.block.length <- max(1, ceiling(ncol(X) / nworkers) * as.numeric(nrow(X)))
      block.length <- min(max.block.length, expected.block.length)
      ## number of columns per block, as calculated by colAutoGrid()
      ncolblock <- min(max(1, floor(block.length / max(1, nrow(X)))), ncol(X))
      grid <- colAutoGrid(X, ncol=.align_to_chunks(ncolblock, X, 2L))
  }
  grid
}

## when the data in 'X' is stored in chunks, such as in HDF5 or Parquet files,
## round down the block width 'n' along the dimension 'whdim' to a multiple of
## the chunk width, so that blocks do not split chunks, which would then be
## read more than once; when a block is narrower than a chunk, its width is
## kept, to avoid exceeding the memory bound behind the block width
#' @importFrom DelayedArray chunkdim
.align_to_chunks <- function(n, X, whdim) {
    cd <- if (is(X, "DelayedArray")) chunkdim(X) else NULL
    if (!is.null(cd) && cd[whdim] > 0 && n >= cd[whdim])
        n <- (n %/% cd[whdim]) * cd[whdim]
    as.integer(n)
}

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitRowsInRanges <- function(grid) {
    rir <- lapply(grid, function(r) ranges(r)[1])
    rir
}

## which elements of the result of bplapply() or bpiterate() did not fail;
## bpiterate() leaves the elements that failed as NULL, and stores their
## errors in the attribute 'errors', named by the position of the element
#' @importFrom BiocParallel bpok
.bp_ok <- function(res) {
    ok <- bpok(res)
    errs <- attr(res, "errors")
    if (!is.null(errs))
        ok[as.integer(names(errs))] <- FALSE

    ok
}

## error of the first element of the result of bplapply() or bpiterate() that
## failed
.bp_first_error <- function(res) {
    errs <- attr(res, "errors")
    if (!is.null(errs))
        return(errs[[1L]])

    res[[which(!bpok(res))[1L]]]
}

#' @importFrom cli cli_alert_warning
.report_parallel_errors <- function(res) {
    bpokmask <- .bp_ok(res)
    msg <- paste("{sum(!bpokmask)} execution thread{?s} gave an error,",
                 "reporting the first one.")
    cli_alert_warning(msg)
    msg <- attr(.bp_first_error(res), "traceback")
    msg <- gsub("\\}", "]", gsub("\\{", "[", msg))
    out <- lapply(msg, cli_alert_warning)
}

## evaluate 'expr', a call to bplapply() or bpiterate(), replacing the
## uninformative error that BiocParallel gives when a forked worker process
## ends without returning its result, such as "wrong args for environment
## subassignment", by an error that tells what happened. that error of
## BiocParallel follows the warning of parallel::mccollect() that a parallel
## job did not deliver its result, whose message is matched in English only, so
## that in other languages the error of BiocParallel is given as it is. errors
## of the workers, which BiocParallel gives as 'bperror' objects, are not
## replaced
#' @importFrom cli cli_abort
.with_dead_worker_check <- function(expr) {
    nodelivery <- FALSE
    withCallingHandlers(expr,
        warning=function(w) {
            if (grepl("parallel jobs? did not deliver", conditionMessage(w)))
                nodelivery <<- TRUE
        },
        error=function(e) {
            if (nodelivery && !inherits(e, "bperror"))
                cli_abort(c("x"=paste("A parallel worker process ended",
                                      "without returning its result."),
                            "i"=paste("This happens, for instance, when the",
                                      "process crashes, e.g., with a",
                                      "segmentation fault, or when the",
                                      "operating system or the workload",
                                      "manager kills it, e.g., for exceeding",
                                      "a memory limit."),
                            "i"=paste("The error output of the process,",
                                      "such as the log of a batch job, may",
                                      "show the cause. Running the same",
                                      "calculations with",
                                      "{.code BPPARAM=SerialParam()} may also",
                                      "show it, but a crash would then end",
                                      "the R session.")),
                          parent=e)
        })
}

## bplapply() and bpiterate() with the arguments in '...', see
## .with_dead_worker_check()
#' @importFrom BiocParallel bplapply
.gsva_bplapply <- function(...) .with_dead_worker_check(bplapply(...))

#' @importFrom BiocParallel bpiterate
.gsva_bpiterate <- function(...) .with_dead_worker_check(bpiterate(...))

## apply 'WRAPPED_FUN' to a block of rows or columns of a matrix sent to a
## worker; defined outside the functions processing matrices by blocks, so
## that sending it to a worker does not send their whole matrix
BLOCK_FUN_WRAPPER <- function(block, WRAPPED_FUN, ...) {
    WRAPPED_FUN(block, ..., verbose=FALSE)
}

## read into main memory the rows (whdim=1) or columns (whdim=2) in the range
## 'rng' of the on-disk matrix 'X', directly through a viewport, because
## subsetting 'X' first, e.g., with X[, rng], takes about a second when 'X'
## binds hundreds of HDF5 datasets, as the output of processing a matrix by
## blocks on disk, see .ondisk_blocks()
#' @importFrom S4Arrays ArrayViewport read_block
#' @importFrom IRanges IRanges start end
.read_block_range <- function(X, whdim, rng) {
    vpstart <- c(1L, 1L)
    vpend <- dim(X)
    vpstart[whdim] <- start(rng)
    vpend[whdim] <- end(rng)
    read_block(X, ArrayViewport(dim(X), IRanges(vpstart, vpend)))
}

## path to a new HDF5 file in the directory 'dumpdir', whose name starts with
## 'prefix' and is unique to the process creating it, see tempfile()
.h5_part_fname <- function(dumpdir, prefix) {
    tempfile(pattern=paste0(prefix, ".", Sys.info()[["nodename"]], ".",
                            Sys.getpid(), "."),
             tmpdir=dumpdir, fileext=".h5")
}

## create in the HDF5 file 'fname', creating it if it does not exist, the
## dataset 'name' of dimensions 'dims', integer or double, and chunks
## 'chunkdim', or no chunks when it has no values. the file is written with
## the rhdf5 package, because the functions of the HDF5Array package that
## create HDF5 datasets lock files shared by all R processes in a way that is
## not safe with concurrent processes
#' @importFrom HDF5Array getHDF5DumpCompressionLevel
#' @importFrom rhdf5 h5createFile h5createDataset
.h5_part_create <- function(fname, dims, int, chunkdim, name="x") {
    if (!file.exists(fname))
        h5createFile(fname)
    empty <- any(dims == 0)
    h5createDataset(fname, name, dims,
                    storage.mode=if (int) "integer" else "double",
                    chunk=if (empty) NULL else pmax(1, pmin(dims, chunkdim)),
                    level=if (empty) 0 else getHDF5DumpCompressionLevel())
}

## number of values of the chunks of the resizable datasets of the CSC layout,
## see .h5_csc_create(). reading a block of columns reads and decompresses
## every chunk holding any of its nonzero values, so that chunks much larger
## than the nonzero values of a block of columns, as in files with a block of
## rows, e.g., the row normalization saved by the jobs of gsvaMap(), make
## reading them by blocks of columns read many times their size, which on
## the shared file systems of computer clusters makes reading the bottleneck
.h5_csc_chunk <- 2^16

## create in the HDF5 file 'fname', creating it if it does not exist, the
## group 'group' holding a sparse
## matrix in compressed sparse column (CSC) layout, as in the HDF5 files of
## 10x Genomics, which the HDF5Array package reads with H5SparseMatrix(): its
## nonzero values "data", integer or double, and their 0-based row indices
## "indices" are resizable datasets, to which the nonzero values of blocks of
## columns are appended, see .h5_csc_append(), while the 0-based positions in
## them of the first nonzero value of each column "indptr" and the dimensions
## "shape" are written at the end, see .h5_csc_close(). it returns the state
## of the appending, see .h5_csc_append()
#' @importFrom HDF5Array getHDF5DumpCompressionLevel
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset H5Sunlimited
.h5_csc_create <- function(fname, int, group="matrix", chunk=.h5_csc_chunk) {
    if (!file.exists(fname))
        h5createFile(fname)
    h5createGroup(fname, group)
    h5createDataset(fname, paste0(group, "/data"), 0, maxdims=H5Sunlimited(),
                    storage.mode=if (int) "integer" else "double",
                    chunk=chunk, level=getHDF5DumpCompressionLevel())
    h5createDataset(fname, paste0(group, "/indices"), 0,
                    maxdims=H5Sunlimited(), storage.mode="integer",
                    chunk=chunk, level=getHDF5DumpCompressionLevel())
    list(fname=fname, group=group, int=int, written=0, nnz=0, indptr=0,
         chunk=chunk, x=if (int) integer(0) else double(0), i=integer(0))
}

## write in the CSC layout the first 'n' nonzero values buffered in 'st', the
## state of the appending, see .h5_csc_append(), returning it updated
#' @importFrom rhdf5 h5set_extent h5write
.h5_csc_flush <- function(st, n) {
    if (n > 0) {
        data <- paste0(st$group, "/data")
        indices <- paste0(st$group, "/indices")
        h5set_extent(st$fname, data, st$written + n)
        h5write(st$x[seq_len(n)], st$fname, data, start=st$written + 1,
                count=n)
        h5set_extent(st$fname, indices, st$written + n)
        h5write(st$i[seq_len(n)], st$fname, indices, start=st$written + 1,
                count=n)
        st$x <- st$x[-seq_len(n)]
        st$i <- st$i[-seq_len(n)]
        st$written <- st$written + n
    }
    st
}

## append to the CSC layout the nonzero values of the sparse block of columns
## 'res', given the state of the appending 'st', see .h5_csc_create(),
## returning it updated. the nonzero values are buffered and written in whole
## chunks, because the HDF5 library does not reuse the space of compressed
## chunks written again, which makes files grow several times their size
.h5_csc_append <- function(st, res) {
    res <- as(res, "CsparseMatrix") ## nonzero values sorted by column
    st$x <- c(st$x, if (st$int) as.integer(res@x) else as.double(res@x))
    st$i <- c(st$i, res@i)
    st$indptr <- c(st$indptr, res@p[-1L] + st$nnz)
    st$nnz <- st$nnz + length(res@i)
    .h5_csc_flush(st, floor(length(st$x) / st$chunk) * st$chunk)
}

## write the nonzero values left in the buffer, the column pointers and the
## dimensions 'dims' of the CSC layout, given the state of the appending 'st',
## see .h5_csc_append()
#' @importFrom rhdf5 h5write
.h5_csc_close <- function(st, dims) {
    st <- .h5_csc_flush(st, length(st$x))
    indptr <- st$indptr
    if (max(indptr) <= .Machine$integer.max)
        indptr <- as.integer(indptr)
    h5write(indptr, st$fname, paste0(st$group, "/indptr"))
    h5write(as.integer(dims), st$fname, paste0(st$group, "/shape"))
}

## apply 'BLOCK_FUN' to a block of rows or columns of a matrix, read into main
## memory when it is stored on disk, so that 'BLOCK_FUN' gives its output in
## main memory; defined outside the functions processing matrices by blocks,
## for the same reason as BLOCK_FUN_WRAPPER()
#' @importFrom S4Arrays DummyArrayGrid read_block
IN_MEMORY_BLOCK_FUN <- function(block, BLOCK_FUN, ..., verbose=FALSE) {
    if (is(block, "DelayedArray"))
        block <- read_block(block, DummyArrayGrid(dim(block))[[1L]])
    BLOCK_FUN(block, ..., verbose=verbose)
}

## apply 'BLOCK_FUN' to a block of rows or columns of a matrix and write its
## output into a new HDF5 file in the directory 'dumpdir', see
## IN_MEMORY_BLOCK_FUN() and .h5_part_fname(), returning the path to that
## file, whether the output is sparse and, if present, its attributes "min"
## and "max". its HDF5 chunks are the default ones of the HDF5Array package
## or, when 'chunkcols' is given, span the rows of the output, with
## 'chunkcols' columns. defined outside the functions processing matrices by
## blocks, for the same reason as BLOCK_FUN_WRAPPER()
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
#' @importFrom HDF5Array getHDF5DumpChunkDim
#' @importFrom rhdf5 h5write
ONDISK_BLOCK_FUN <- function(block, BLOCK_FUN, dumpdir, prefix,
                             chunkcols=NULL, ..., verbose=FALSE) {
    res <- IN_MEMORY_BLOCK_FUN(block, BLOCK_FUN, ..., verbose=verbose)
    rm(block)

    fname <- .h5_part_fname(dumpdir, prefix)
    sparse <- is_sparse(res)
    int <- type(res) == "integer"
    rmin <- attr(res, "min")
    rmax <- attr(res, "max")
    res <- as.matrix(res) ## HDF5 datasets are dense
    attributes(res) <- list(dim=dim(res))
    chunkdim <- getHDF5DumpChunkDim(dim(res))
    if (!is.null(chunkcols))
        chunkdim <- c(nrow(res), chunkcols)
    .h5_part_create(fname, dim(res), int, chunkdim)
    h5write(res, fname, "x")

    list(fname=fname, sparse=sparse, min=rmin, max=rmax)
}

## apply 'BLOCK_FUN' to each of the consecutive blocks of rows (whdim=1) or
## columns (whdim=2) of a matrix with the ranges in 'grp', taken from the
## matrix and processed in main memory by 'FUN_WRAPPER', see .bp_blocks()
## and IN_MEMORY_BLOCK_FUN(), and write their outputs into a single new HDF5
## file in the directory 'dumpdir', see .h5_part_fname(), returning the path
## to that file, whether the output is sparse and stored in CSC layout and, if
## present, the minimum and maximum of its attributes "min" and "max". sparse
## output of blocks of columns, such as sparse ranks, is stored in CSC layout,
## see .h5_csc_create(), which takes less space and is read faster by blocks
## of columns than a dense HDF5 dataset, as measured on single-cell data.
## otherwise, the output is stored in a dense HDF5 dataset whose chunks span
## the blocks, so that each block fills whole chunks: when the blocks are of
## rows, chunks span their rows, with 'chunkcols' columns, and when they are
## of columns, chunks span their columns, with as many rows as the default
## chunk length of the HDF5Array package takes. 'verbose' and 'idpbe' are given to
## 'FUN_WRAPPER' to report progress. defined outside the functions processing
## matrices by blocks, for the same reason as BLOCK_FUN_WRAPPER()
#' @importFrom S4Arrays is_sparse
#' @importFrom BiocGenerics type
#' @importFrom HDF5Array getHDF5DumpChunkLength
#' @importFrom rhdf5 h5write
#' @importFrom IRanges width
ONDISK_GROUP_FUN <- function(grp, FUN_WRAPPER, BLOCK_FUN, whdim, dumpdir,
                             prefix, chunkcols=NULL, ..., verbose=FALSE,
                             idpbe=NULL) {
    fname <- .h5_part_fname(dumpdir, prefix)
    total <- sum(vapply(grp, width, integer(1)))
    offset <- 0L
    sparse <- csc <- NULL
    mines <- maxes <- NULL
    for (rng in grp) {
        ## the memory of the previous block, its copy taken from the input
        ## and its output, is released before processing the next one,
        ## because otherwise R collects it only after the memory allocated
        ## reaches a threshold, which grows with the memory used by the
        ## input when it is in main memory, piling up several blocks
        res <- NULL
        invisible(gc(full=FALSE))
        res <- FUN_WRAPPER(rng, verbose=verbose, idpbe=idpbe,
                           WRAPPED_FUN=IN_MEMORY_BLOCK_FUN, BLOCK_FUN=BLOCK_FUN,
                           ...)
        if (!is.null(attr(res, "min"))) {
            mines <- min(mines, attr(res, "min"))
            maxes <- max(maxes, attr(res, "max"))
        }
        if (is.null(sparse)) { ## the file is created with the first block
            sparse <- is_sparse(res)
            csc <- sparse && whdim == 2L
            int <- type(res) == "integer"
            dims <- dim(res)
            dims[whdim] <- total
            if (csc)
                st <- .h5_csc_create(fname, int)
            else {
                if (whdim == 1L)
                    chunkdim <- c(nrow(res), chunkcols)
                else
                    chunkdim <- c(floor(getHDF5DumpChunkLength() / ncol(res)),
                                  ncol(res))
                .h5_part_create(fname, dims, int, chunkdim)
            }
        }
        if (csc) {
            st <- .h5_csc_append(st, res)
            next
        }
        res <- as.matrix(res) ## HDF5 datasets are dense
        attributes(res) <- list(dim=dim(res))
        start <- c(1L, 1L)
        start[whdim] <- offset + 1L
        h5write(res, fname, "x", start=start, count=dim(res))
        offset <- offset + dim(res)[whdim]
    }
    if (csc)
        .h5_csc_close(st, dims)

    list(fname=fname, sparse=sparse, csc=csc, min=mines, max=maxes)
}

## split the list of consecutive ranges 'rngs' into at most 'ngroups' groups
## of consecutive ranges with about the same number of ranges
.group_ranges <- function(rngs, ngroups) {
    ngroups <- max(1L, min(length(rngs), ngroups))
    unname(split(rngs, ceiling(seq_along(rngs) * ngroups / length(rngs))))
}

## process in blocks of rows (whdim=1) or columns (whdim=2), with the ranges
## 'rngs', the matrix 'X' with 'FUN', giving its output in an on-disk data
## structure, see .processMatrixRows() and .processMatrixCols(). the blocks
## are split into one group of consecutive blocks per worker, whose outputs
## are written into one HDF5 file per group, see ONDISK_GROUP_FUN(), and the
## output is formed by binding the HDF5 data of those files, because reading
## blocks out of many small files is slower, as measured on single-cell data,
## where more groups per worker were not faster; serially, there is a single
## group and a single file. workers
## not forked from this process processing a matrix in main memory receive
## its blocks one at a time, see .bp_blocks(), so that the output of each
## block is written into a separate HDF5 file, see ONDISK_BLOCK_FUN(). when
## the blocks are of rows, the chunks of those files are as wide as the
## default blocks of columns of the output, because it is later read by
## blocks of columns, e.g., to rank its columns, and wider chunks would be
## read and decompressed many times. 'FUN_WRAPPER' takes the blocks from 'X'
## in the workers. the minimum and maximum enrichment scores of ssGSEA,
## stored in the attributes "min" and "max" of the output of each block, are
## kept in the output
#' @importFrom cli cli_abort cli_alert_warning cli_progress_bar
#' @importFrom cli cli_progress_done
#' @importFrom BiocParallel bptry "bpprogressbar<-" bpprogressbar
#' @importFrom BiocParallel bpstopOnError "bpstopOnError<-" MulticoreParam
#' @importFrom HDF5Array HDF5Array H5SparseMatrix getHDF5DumpDir
#' @importFrom DelayedArray getAutoBlockLength
.ondisk_blocks <- function(X, whdim, rngs, FUN, FUN_WRAPPER, ..., nworkers,
                           verbose, progressmsg, BPPARAM) {
    dumpdir <- getHDF5DumpDir()
    if (!dir.exists(dumpdir))
        dir.create(dumpdir, recursive=TRUE)
    prefix <- basename(tempfile("GSVA")) ## unique to this call
    chunkcols <- NULL
    if (whdim == 1L) ## the output is double, as the one of normalizing rows
        chunkcols <- max(1, floor(getAutoBlockLength("double") /
                                  max(1, nrow(X))))

    if (nworkers <= 1L) {
        env <- NULL
        if (verbose) {
            env <- new.env(parent=globalenv())
            assign("idpb", cli_progress_bar(progressmsg, total=dim(X)[whdim]),
                   envir=env)
        }
        res <- list(ONDISK_GROUP_FUN(rngs, FUN_WRAPPER, FUN, whdim, dumpdir,
                                     prefix, chunkcols, ..., verbose=verbose,
                                     idpbe=env))
        if (verbose)
            cli_progress_done(get("idpb", envir=env))
    } else {
        ## 'BPPARAM' is a reference class object, so the following changes
        ## reach the object of the caller, which is restored on exit
        oldprogressbar <- bpprogressbar(BPPARAM)
        oldstoponerror <- bpstopOnError(BPPARAM)
        on.exit({
            bpprogressbar(BPPARAM) <- oldprogressbar
            bpstopOnError(BPPARAM) <- oldstoponerror
        }, add=TRUE)
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bpstopOnError(BPPARAM) <- FALSE
        grouped <- is(X, "DelayedArray") || is(BPPARAM, "MulticoreParam")
        run <- function(BPREDO=list()) {
            if (grouped)
                .gsva_bplapply(.group_ranges(rngs, nworkers),
                               FUN=ONDISK_GROUP_FUN, FUN_WRAPPER=FUN_WRAPPER,
                               BLOCK_FUN=FUN, whdim=whdim, dumpdir=dumpdir,
                               prefix=prefix, chunkcols=chunkcols, ...,
                               BPREDO=BPREDO, BPPARAM=BPPARAM)
            else
                .bp_blocks(X, whdim, rngs, ONDISK_BLOCK_FUN, FUN_WRAPPER,
                           BLOCK_FUN=FUN, dumpdir=dumpdir, prefix=prefix,
                           chunkcols=chunkcols, ..., BPREDO=BPREDO,
                           BPPARAM=BPPARAM)
        }
        res <- bptry(run())
        if (any(!.bp_ok(res))) {
            .report_parallel_errors(res)
            cli_alert_warning("Trying to execute again the failing thread(s)")
            res <- bptry(run(BPREDO=res))
            if (any(!.bp_ok(res))) {
                .report_parallel_errors(res)
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }

    pieces <- lapply(res, function(r)
        if (isTRUE(r$csc)) H5SparseMatrix(r$fname, "matrix") else
        HDF5Array(r$fname, "x", as.sparse=r$sparse))
    out <- if (length(pieces) == 1L) pieces[[1L]] else
           do.call(if (whdim == 1L) "rbind" else "cbind", pieces)

    mines <- unlist(lapply(res, "[[", "min"))
    if (!is.null(mines)) { ## min and max enrichment scores stored by ssGSEA
        attr(out, "min") <- min(mines)
        attr(out, "max") <- max(unlist(lapply(res, "[[", "max")))
    }

    out
}

## blocks of rows (whdim=1) or columns (whdim=2) of the matrix 'X' in main
## memory, with the ranges 'rngs', processed by the workers of 'BPPARAM'. workers
## not forked from this process, such as socket workers, receive the blocks
## one at a time, as they become free, while forked workers, and workers
## processing a matrix stored on disk, take their blocks from 'X' themselves
## with 'FUN_WRAPPER', which would send the whole 'X' to workers not forked
## from this process. 'BPREDO' gives the result of a previous call to
## recompute only its failed elements
#' @importFrom BiocParallel MulticoreParam
#' @importFrom IRanges start end
.bp_blocks <- function(X, whdim, rngs, FUN, FUN_WRAPPER, ..., BPREDO=list(),
                       BPPARAM) {
    if (is(X, "DelayedArray") || is(BPPARAM, "MulticoreParam"))
        return(.gsva_bplapply(rngs, FUN=FUN_WRAPPER, verbose=FALSE,
                              idpbe=NULL, WRAPPED_FUN=FUN, ..., BPREDO=BPREDO,
                              BPPARAM=BPPARAM))

    i <- 0L
    ITER <- function() {
        i <<- i + 1L
        if (i > length(rngs))
            return(NULL)
        rng <- start(rngs[[i]]):end(rngs[[i]])
        if (whdim == 1L) X[rng, , drop=FALSE] else X[, rng, drop=FALSE]
    }
    .gsva_bpiterate(ITER, FUN=BLOCK_FUN_WRAPPER, WRAPPED_FUN=FUN, ...,
                    BPREDO=BPREDO, BPPARAM=BPPARAM)
}

## replace in the seeds of the on-disk matrix 'X' the H5File objects of the
## h5mread package, which hold an open connection to an HDF5 file, by the
## path to that file, when the file is local. the connection of an H5File
## object is lost when it is serialized, e.g., to send it to the workers of a
## SnowParam back-end, and it is not reliable in the workers of a
## MulticoreParam back-end, as the documentation of H5File explains, while a
## file path is opened by each process reading from it. a remote file, such
## as one read from Amazon S3, can only be read through its H5File object, so
## that processing it with 'nworkers' workers of a back-end that does not
## fork this process, such as SnowParam, gives an error before starting
#' @importFrom DelayedArray modify_seeds seedApply
#' @importFrom cli cli_abort
.h5file_seeds_to_paths <- function(X, BPPARAM, nworkers) {
    if (!is(X, "DelayedArray"))
        return(X)
    h5file <- function(s) .hasSlot(s, "filepath") && is(s@filepath, "H5File")
    if (!any(unlist(seedApply(X, h5file))))
        return(X)

    X <- modify_seeds(X, function(s) {
        if (h5file(s) && !isTRUE(s@filepath@s3))
            s@filepath <- s@filepath@filepath
        s
    })
    if (nworkers > 1L && !is(BPPARAM, "MulticoreParam") &&
        any(unlist(seedApply(X, h5file)))) {
        bpclass <- class(BPPARAM)[1]
        cli_abort(c("x"=paste("The input data is read from a remote HDF5",
                              "file through an {.cls H5File} object, which",
                              "cannot be sent to the parallel workers of a",
                              "{.cls {bpclass}} back-end, see",
                              "{.help h5mread::H5File}."),
                    "i"=paste("Use {.code BPPARAM=SerialParam()}, or a",
                              "local copy of the HDF5 file.")))
    }

    X
}

## process the rows of a matrix with a given function FUN, opening parallelism
## through a BiocParallelParam object BPPARAM, when different from NULL, and
## reporting progress using the 'cli' package when possible. the arguments
## 'workfactor', 'outfactor' and 'outextra' give the memory that FUN takes to
## process the rows of 'X' within the maximum main memory 'maxmem', see
## .units_per_block(), as do 'dense', when FUN converts blocks of sparse 'X'
## into dense matrices, and 'heldmem', the memory held by other objects
## with 'sinkout=TRUE', FUN processes each block in main memory and the
## output is written into on-disk data structures, see .ondisk_blocks()

#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort cli_progress_bar cli_alert_warning
#' @importFrom BiocParallel bplapply bpnworkers "bpprogressbar<-" bptry
#' @importFrom BiocParallel bpok bpstopOnError "bpstopOnError<-"
#' @importFrom memuse howbig
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixRows <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf, workfactor=2, outfactor=1,
                               outextra=0, dense=FALSE,
                               heldmem=.held_mem(), sinkout=FALSE) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        stopifnot(is(BPPARAM, "BiocParallelParam"))
        nworkers <- bpnworkers(BPPARAM)
    }
    X <- .h5file_seeds_to_paths(X, BPPARAM, nworkers)

    grid <- .rowgridsize(X, nworkers, maxmem, workfactor, outfactor, outextra,
                         dense=dense, heldmem=heldmem)
    rir <- .splitRowsInRanges(grid)
    if (length(rir) > 1 && verbose) {
        sze <- howbig(as.numeric(width(rir[[1]])), as.numeric(ncol(X)),
                      representation="dense", type=type(X))
        msg <- sprintf("Splitting calculations in %d chunks of [%d, %d] and %s",
                       length(rir), width(rir[[1]]), ncol(X), as.character(sze))
        cli_alert_info(msg)
    } else if (length(rir) == 1 && !sinkout) ## serial execution in one call
        return(FUN(X, ..., verbose=verbose))

    FUN_WRAPPER <- function(rowsrng, verbose, idpbe, WRAPPED_FUN, ...) {
        rng <- rowsrng
        if (!is(X, "DelayedMatrix"))
            rng <- start(rowsrng):end(rowsrng)
        block <- if (sinkout && is(X, "DelayedMatrix"))
                     .read_block_range(X, 1L, rowsrng)
                 else
                     X[rng, , drop=FALSE]
        res <- WRAPPED_FUN(block, ..., verbose=FALSE)
        if (verbose && is(idpbe, "environment"))
            cli_progress_update(id=get("idpb", envir=idpbe), width(rowsrng))
        return(res)
    }

    if (sinkout)
        return(.ondisk_blocks(X, 1L, rir, FUN, FUN_WRAPPER, ...,
                              nworkers=nworkers, verbose=verbose,
                              progressmsg=progressmsg, BPPARAM=BPPARAM))
        
    totalnrows <- nrow(X)
    res <- NULL
    if (is.null(BPPARAM) || nworkers <= 1L) { ## serial execution in chunks
        env <- NULL
        if (verbose) {
            env <- new.env(parent=globalenv())
            assign("idpb", cli_progress_bar(progressmsg, total=totalnrows),
                   envir=env)
        }
        res <- lapply(rir, FUN=FUN_WRAPPER, verbose=verbose,
                      idpbe=env, WRAPPED_FUN=FUN, ...)
        if (verbose)
            cli_progress_done(get("idpb", envir=env))
    } else {                                  ## parallel execution in chunks
        ## 'BPPARAM' is a reference class object, so the following changes
        ## reach the object of the caller, which is restored on exit
        oldprogressbar <- bpprogressbar(BPPARAM)
        oldstoponerror <- bpstopOnError(BPPARAM)
        on.exit({
            bpprogressbar(BPPARAM) <- oldprogressbar
            bpstopOnError(BPPARAM) <- oldstoponerror
        }, add=TRUE)
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bpstopOnError(BPPARAM) <- FALSE
        res <- bptry(.bp_blocks(X, 1L, rir, FUN, FUN_WRAPPER, ...,
                                BPPARAM=BPPARAM))
        bpokmask <- .bp_ok(res)
        if (any(!bpokmask)) {
            .report_parallel_errors(res)
            cli_alert_warning("Trying to execute again the failing thread(s)")
            res <- bptry(.bp_blocks(X, 1L, rir, FUN, FUN_WRAPPER, ...,
                                    BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- .bp_ok(res)
            if (any(!bpokmask)) {
                .report_parallel_errors(res)
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }
    res <- do.call("rbind", res)

    return(res)
}

#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges ranges
.splitColsInRanges <- function(grid) {
    cir <- lapply(grid, function(r) ranges(r)[2])
    cir
}

## bind the blocks of columns in the list 'blocks'. dgCMatrix blocks are bound
## by concatenating their slots, which allocates only the result, while
## binding them with cbind() allocates intermediate matrices that take, as
## measured, from two to four times the size of the result as the number of
## blocks grows from 4 to 40
.cbind_blocks <- function(blocks) {
    if (length(blocks) < 2L ||
        !all(vapply(blocks, is, logical(1), "dgCMatrix")))
        return(do.call("cbind", blocks))

    nzpercol <- unlist(lapply(blocks, function(b) diff(b@p)), use.names=FALSE)
    colnms <- unlist(lapply(blocks, colnames), use.names=FALSE)
    new("dgCMatrix",
        Dim=c(nrow(blocks[[1L]]), length(nzpercol)),
        Dimnames=list(rownames(blocks[[1L]]), colnms),
        i=unlist(lapply(blocks, slot, "i"), use.names=FALSE),
        x=unlist(lapply(blocks, slot, "x"), use.names=FALSE),
        p=c(0L, cumsum(nzpercol)))
}

## process the columns of a matrix with a given function FUN, opening parallelism
## through a BiocParallelParam object BPPARAM, when different from NULL, and
## reporting progress using the 'cli' package when possible. the arguments
## 'workfactor', 'outfactor' and 'outextra' give the memory that FUN takes to
## process the columns of 'X' within the maximum main memory 'maxmem', see
## .units_per_block(), as do 'dense', when FUN converts blocks of sparse 'X'
## into dense matrices, and 'heldmem', the memory held by other objects
## the results of FUN on each block of columns are bound by columns or,
## when 'combine' is a function of two results, reduced with it
## with 'sinkout=TRUE', FUN processes each block in main memory and the
## output is written into on-disk data structures, see .ondisk_blocks()
#' @importFrom BiocGenerics type
#' @importFrom cli cli_abort
#' @importFrom BiocParallel bplapply bpnworkers "bpprogressbar<-" bptry
#' @importFrom BiocParallel bpok bpstopOnError "bpstopOnError<-"
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges start end width
.processMatrixCols <- function(X, FUN, ..., verbose=TRUE,
                               minparrows=100, minparcols=100,
                               progressmsg="Progress", BPPARAM=NULL,
                               maxmem=Inf, workfactor=2, outfactor=1,
                               outextra=0, dense=FALSE,
                               heldmem=.held_mem(), combine=NULL,
                               sinkout=FALSE) {
    stopifnot(length(dim(X)) == 2) ## QC
    FUN <- match.fun(FUN)
    nworkers <- 1L
    if (!is.null(BPPARAM) && nrow(X) > minparrows && ncol(X) > minparcols) {
        stopifnot(is(BPPARAM, "BiocParallelParam"))
        nworkers <- bpnworkers(BPPARAM)
    }
    X <- .h5file_seeds_to_paths(X, BPPARAM, nworkers)

    grid <- .colgridsize(X, nworkers, maxmem, workfactor, outfactor, outextra,
                         dense=dense, heldmem=heldmem)
    cir <- .splitColsInRanges(grid)

    if (length(cir) > 1 && verbose) {
        sze <- howbig(as.numeric(nrow(X)), as.numeric(width(cir[[1]])),
                      representation="dense", type=type(X))
        cli_alert_info(sprintf("Splitting calculations in %d chunks of [%d, %d] and %s",
                               length(cir), nrow(X), width(cir[[1]]), as.character(sze)))
    } else if (length(cir) == 1 && !sinkout) ## serial execution in one call
        return(FUN(X, ..., verbose=verbose))

    FUN_WRAPPER <- function(colsrng, verbose, idpbe, WRAPPED_FUN, ...) {
        rng <- colsrng
        if (!is(X, "DelayedMatrix"))
            rng <- start(colsrng):end(colsrng)
        block <- if (sinkout && is(X, "DelayedMatrix"))
                     .read_block_range(X, 2L, colsrng)
                 else
                     X[, rng, drop=FALSE]
        res <- WRAPPED_FUN(block, ..., verbose=FALSE)
        if (verbose && is(idpbe, "environment"))
            cli_progress_update(id=get("idpb", envir=idpbe), width(colsrng))
        return(res)
    }

    if (sinkout)
        return(.ondisk_blocks(X, 2L, cir, FUN, FUN_WRAPPER, ...,
                              nworkers=nworkers, verbose=verbose,
                              progressmsg=progressmsg, BPPARAM=BPPARAM))
        
    totalncols <- ncol(X)
    res <- NULL
    if (is.null(BPPARAM) || nworkers <= 1L) { ## serial execution in chunks
        env <- new.env(parent=globalenv())
        assign("idpb", cli_progress_bar(progressmsg, total=totalncols), envir=env)
        res <- lapply(cir, FUN=FUN_WRAPPER, verbose=verbose,
                      idpbe=env, WRAPPED_FUN=FUN, ...)
        cli_progress_done(get("idpb", envir=env))
    } else {                                  ## parallel execution in chunks
        ## 'BPPARAM' is a reference class object, so the following changes
        ## reach the object of the caller, which is restored on exit
        oldprogressbar <- bpprogressbar(BPPARAM)
        oldstoponerror <- bpstopOnError(BPPARAM)
        on.exit({
            bpprogressbar(BPPARAM) <- oldprogressbar
            bpstopOnError(BPPARAM) <- oldstoponerror
        }, add=TRUE)
        if (verbose)
            bpprogressbar(BPPARAM) <- TRUE    ## reporting progress wo/ cli
        bpstopOnError(BPPARAM) <- FALSE
        res <- bptry(.bp_blocks(X, 2L, cir, FUN, FUN_WRAPPER, ...,
                                BPPARAM=BPPARAM))
        bpokmask <- .bp_ok(res)
        if (any(!bpokmask)) {
            .report_parallel_errors(res)
            cli_alert_warning("Trying to execute again the failing thread(s)")
            res <- bptry(.bp_blocks(X, 2L, cir, FUN, FUN_WRAPPER, ...,
                                    BPREDO=res, BPPARAM=BPPARAM))
            bpokmask <- .bp_ok(res)
            if (any(!bpokmask)) {
                .report_parallel_errors(res)
                cli_abort(c("x"="Cancelling execution"))
            }
        }
    }

    if (is.function(combine))
        return(Reduce(combine, res))

    mines <- maxes <- NULL
    if (!is.null(attr(res[[1]], "min"))) { ## min and max enrichment scores stored by ssGSEA
        mines <- min(vapply(X=res, FUN=function(x) attr(x, "min"), FUN.VALUE=numeric(1)))
        maxes <- max(vapply(X=res, FUN=function(x) attr(x, "max"), FUN.VALUE=numeric(1)))
    }

    res <- .cbind_blocks(res)

    if (!is.null(mines)) { ## min and max enrichment scores stored by ssGSEA
        attr(res, "min") <- mines
        attr(res, "max") <- maxes
    }

    return(res)
}
