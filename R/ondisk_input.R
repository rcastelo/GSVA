## This file contains the internal functions that load into main memory sparse
## input data stored on disk, by blocks of columns read in parallel when
## possible, and that count or estimate the number of nonzero values of the
## input data, from which its size in main memory is estimated. They are not
## exported and are only used internally by GSVA.

## load into main memory as an SVT_SparseMatrix object the sparse on-disk data
## in 'X', reading one block of columns at a time, because coercing 'X' as a
## whole with as() takes 5 to 8 times the memory of the result. reading a
## block as sparse with the HDF5Array package takes memory while it lasts,
## besides the blocks of columns already read, which form the result and take
## their memory when bound: when 'X' is stored in chunks of several columns,
## as in dense HDF5 datasets, each chunk is decompressed into its dense form
## before its nonzero values are taken, so that a block takes a multiple of
## the size of its dense form, whatever its fraction of nonzero values, see
## 'densereadfactor' in .memory_model_defaults, and otherwise, as in the CSC
## layout, a number of bytes per nonzero value, see 'sparsereadbytes'. when
## 'X' is stored in chunks of several columns, also through delayed
## operations such as a subset of its columns, see .ondisk_chunkdim(), blocks
## span the width of the chunks,
## because blocks splitting chunks make their reading slower, while wider
## blocks do not make it faster: the smallest multiple of that width not
## narrower than the default block of the DelayedArray package. otherwise, as
## in the compressed
## sparse column (CSC) layout of HDF5 files, e.g., those written by
## saveHDF5GSVA(), where each column is a chunk, or when 'X' binds several
## files, blocks are as wide as the memory allows, because each block reads
## each file anew, which takes a long time on the shared file systems of
## computer clusters, as measured with the row normalization saved by the
## jobs of gsvaMap(); without knowing that memory, blocks have the default
## size. blocks are narrower when their memory does not fit in the fraction
## .mem_fraction_R() of the maximum main memory 'maxmem' left by the result,
## whose estimated size is 'insize', when both are known, by the memory held
## by other objects, see .held_mem(), and by the memory of this R process by
## itself, see .block_worker_mem(), but not narrower than a chunk, or
## than the default block when 'X' is not stored in chunks of several columns,
## because narrower blocks would make the loading take much longer, e.g.,
## blocks of a single column would decompress all the chunks that it spans for
## each column, when the memory left by the result is too small, where the
## whole calculations do not fit in 'maxmem' anyway. the number of nonzero
## values of
## each column, when the layout of 'X' gives it without reading its values,
## see .col_nzcounts(), sets the width of each block to the columns that fit
## in that memory, because columns with many more nonzero values than the
## average, e.g., of cells with deeper sequencing, may be stored next to each
## other, and the average number of nonzero values per column, estimated
## from 'insize', otherwise. with the parallel workers of 'BPPARAM', the
## blocks are read by the workers, see .load_workers(), when 'X' is stored in
## chunks of several columns, while otherwise they are read in this process,
## because reading the CSC layout of 100,000 cells of single-cell data took
## about 20% to 45% longer with six workers than in this process, on the
## shared file system of a computer cluster, and twice as long on macOS,
## since each column is read without decompressing other columns and sending
## the blocks back from the workers takes longer than reading them; with
## 'parallelonly=TRUE', it gives NULL instead of reading them in this process
## when fewer than two workers can read them, see .load_dense()
#' @importFrom DelayedArray colAutoGrid read_block ArbitraryArrayGrid
#' @importFrom BiocGenerics type
.load_sparse_by_blocks <- function(X, maxmem=Inf, insize=NA_real_,
                                   BPPARAM=NULL, parallelonly=FALSE) {
    if (ncol(X) == 0L)
        return(as(X, "SVT_SparseMatrix"))

    ncolblock <- ncol(colAutoGrid(X)[[1L]]) ## default block size
    budget <- !is.null(maxmem) && !is.null(insize) && is.finite(maxmem) &&
              !is.na(insize)
    cd <- .ondisk_chunkdim(X)
    multicol <- !is.null(cd) && cd[2] > 1
    if (multicol) ## multiple of the width of the chunks
        ncolblock <- cd[2] * ceiling(ncolblock / cd[2])
    ## the narrowest blocks within the memory available
    mincols <- if (multicol) cd[2] else ncolblock
    if (!multicol && budget) ## as wide as the memory allows, see below
        ncolblock <- ncol(X)
    ## the DelayedArray package does not support blocks with more than
    ## .Machine$integer.max values
    maxcols <- floor(.Machine$integer.max / max(1, nrow(X)))
    grid <- NULL
    avail <- Inf
    blockbytes <- NA_real_ ## memory of the largest block while it is read
    pm <- NULL
    nzbytes <- .memory_model("sparsereadbytes")
    eltbytes <- if (type(X) == "integer") 4 else 8
    if (budget) {
        pm <- .block_worker_mem(X, 2L, forked=.forked_workers(BPPARAM))
        avail <- .mem_fraction_R() * maxmem - insize - .held_mem() - pm$workermem
        nz <- if (multicol) NULL else .col_nzcounts(X)
        if (!is.null(nz)) {
            ## at least the memory of 'mincols' columns of average width
            avail <- max(avail, nzbytes * mean(pmax(1, nz)) * mincols)
            grid <- .nzcount_col_grid(X, nzbytes * pmax(1, nz), avail,
                                      min(ncolblock, maxcols))
            ends <- cumsum(vapply(seq_along(grid), function(i)
                                      ncol(grid[[i]]), integer(1)))
            blockbytes <- max(diff(c(0, cumsum(nzbytes *
                                               pmax(1, nz))[ends])))
        } else {
            ## memory of each column while it is read, by the size of its
            ## dense form when it is stored in chunks of several columns, and
            ## otherwise by its average number of nonzero values
            nzpercol <- insize / (eltbytes + 4) / ncol(X)
            colbytes <- if (multicol) .memory_model("densereadfactor") *
                                      as.numeric(nrow(X)) * eltbytes
                        else nzbytes * max(1, nzpercol)
            ncolblock <- max(mincols,
                             min(ncolblock, floor(avail / colbytes)))
            if (multicol) ## whole chunks
                ncolblock <- cd[2] * (ncolblock %/% cd[2])
        }
    }
    if (is.null(grid)) {
        ncolblock <- as.integer(max(1, min(ncol(X), ncolblock, maxcols)))
        grid <- colAutoGrid(X, ncol=ncolblock)
        if (budget)
            blockbytes <- colbytes * ncolblock
    }
    nworkers <- if (multicol) .load_workers(X, length(grid), avail,
                                            blockbytes, BPPARAM, pm)
                else 1L
    if (parallelonly && nworkers < 2L)
        return(NULL)
    if (nworkers > 1L)
        blocks <- .read_sparse_blocks_bp(X, grid, nworkers, BPPARAM)
    else
        blocks <- lapply(seq_along(grid), .read_sparse_block, X, grid)
    res <- do.call(cbind, blocks)
    dimnames(res) <- dimnames(X)

    res
}

## block 'i' of the grid 'grid' of the on-disk sparse matrix 'M', read into
## main memory as an SVT_SparseMatrix object, see .load_sparse_by_blocks();
## the matrix is not called 'X', which is the first argument of bplapply()
#' @importFrom DelayedArray read_block
.read_sparse_block <- function(i, M, grid) {
    as(read_block(M, grid[[i]], as.sparse=TRUE), "SVT_SparseMatrix")
}

## number of parallel workers of 'BPPARAM' that read at the same time the
## 'nblocks' blocks of the on-disk sparse matrix 'X', each one taking at most
## 'blockbytes' bytes while it is read, within the memory 'avail' left by the
## result, see .load_sparse_by_blocks(): as many as there are workers and
## blocks, when the memory is not known, and otherwise as many as fit in
## 'avail' together with the memory of all the workers, or 1, to read the
## blocks in this process, when fewer than two fit. forked workers, such as
## those of a MulticoreParam object, return each block as soon as it is read,
## without keeping memory of their own from one block to the next, so that
## each one takes, besides its block, the memory outside R of a forked
## worker with its input in main memory, see .forked_worker_mem(), and its
## buffers to read 'X' from disk, see .disk_read_mem(): together with a block
## of 324 MiB in its dense form, at most about 0.75 GiB on Linux and 1.5 GiB
## on macOS, as measured with two to six workers reading 100,000 to 300,000
## cells;
## workers not forked take the memory of .block_worker_mem(), given in 'pm'.
## reading the blocks of 300,000 cells of single-cell data stored in a dense
## HDF5 dataset took a third of the time with six workers
#' @importFrom BiocParallel bpnworkers
.load_workers <- function(X, nblocks, avail, blockbytes, BPPARAM,
                          pm=NULL) {
    if (is.null(BPPARAM) || nblocks < 2L)
        return(1L)
    nworkers <- min(bpnworkers(BPPARAM), nblocks)
    if (nworkers < 2L || !is.finite(avail) || is.na(blockbytes))
        return(as.integer(nworkers))

    if (.forked_workers(BPPARAM))
        procmem <- .forked_worker_mem(inmemory=TRUE) + .disk_read_mem()
    else {
        if (is.null(pm))
            pm <- .block_worker_mem(X, 2L, forked=FALSE)
        procmem <- pm$parworkermem + pm$copiedmem
    }
    procs <- bpnworkers(BPPARAM) * procmem
    as.integer(max(1, min(nworkers, floor((avail - procs) / blockbytes))))
}

## blocks of the grid 'grid' of the on-disk sparse matrix 'X' read by the
## workers of 'BPPARAM', in rounds of 'nworkers' blocks, so that at most
## 'nworkers' blocks are read at the same time, and each worker returns each
## of its blocks as soon as it is read, instead of holding all of them until
## the end, as bplapply() does with the elements given to a worker. the
## workers run all the rounds, instead of being started anew for each of
## them, and without progress bars, see .bp_blocks()
#' @importFrom BiocParallel bpisup bpstart bpstop bpprogressbar
#' @importFrom BiocParallel "bpprogressbar<-"
.read_sparse_blocks_bp <- function(X, grid, nworkers, BPPARAM) {
    X <- .h5file_seeds_to_paths(X, BPPARAM, nworkers)
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }
    if (bpprogressbar(BPPARAM)) {
        bpprogressbar(BPPARAM) <- FALSE
        on.exit(bpprogressbar(BPPARAM) <- TRUE, add=TRUE)
    }
    rounds <- split(seq_along(grid),
                    ceiling(seq_along(grid) / nworkers))
    unlist(lapply(rounds, function(r)
                      .gsva_bplapply(r, .read_sparse_block, M=X, grid=grid,
                                     BPPARAM=BPPARAM)),
           recursive=FALSE, use.names=FALSE)
}

## number of nonzero values of each column of the sparse matrix 'X', when
## they can be known without reading its values: from the pointers to the
## start of each column in the compressed sparse column (CSC) layout of
## HDF5 files, e.g., those written by saveHDF5GSVA() or by 10x Genomics,
## through the delayed operations that restrict, bind along the columns or
## name them; those of whole columns bound those of a subset of their rows.
## NULL otherwise
#' @importClassesFrom HDF5Array CSC_H5SparseMatrixSeed
#' @importClassesFrom DelayedArray DelayedSubset DelayedAbind DelayedSetDimnames
.col_nzcounts <- function(X) {
    nz <- .seed_col_nzcounts(if (is(X, "DelayedArray")) X@seed else X)
    if (is.null(nz) || length(nz) != ncol(X))
        return(NULL)

    as.numeric(nz)
}

.seed_col_nzcounts <- function(s) {
    if (is(s, "CSC_H5SparseMatrixSeed"))
        return(s@indptr_ranges$width)
    if (is(s, "DelayedArray") || is(s, "DelayedSetDimnames"))
        return(.seed_col_nzcounts(s@seed))
    if (is(s, "DelayedSubset")) {
        nz <- .seed_col_nzcounts(s@seed)
        j <- s@index[[2L]]
        if (is.null(nz) || is.null(j))
            return(nz)
        return(nz[j])
    }
    if (is(s, "DelayedAbind") && s@along == 2L) {
        nzs <- lapply(s@seeds, .seed_col_nzcounts)
        if (any(vapply(nzs, is.null, logical(1))))
            return(NULL)
        return(unlist(nzs, use.names=FALSE))
    }

    NULL
}

## grid of blocks of consecutive columns of 'X', whose columns take 'bytes'
## bytes each while read, of at most 'maxcols' columns and as many columns as
## fit in 'avail' bytes, or a single column when it does not fit by itself
#' @importFrom DelayedArray ArbitraryArrayGrid
.nzcount_col_grid <- function(X, bytes, avail, maxcols) {
    n <- length(bytes)
    cs <- cumsum(bytes)
    ends <- integer(0)
    last <- 0L
    while (last < n) {
        used <- if (last == 0L) 0 else cs[last]
        end <- findInterval(used + avail, cs)
        end <- min(max(end, last + 1L), last + maxcols, n)
        ends <- c(ends, as.integer(end))
        last <- end
    }

    ArbitraryArrayGrid(list(nrow(X), ends))
}

## load into main memory as a dense matrix the on-disk data in 'X'. with the
## parallel workers of 'BPPARAM', when 'X' is stored in chunks of several
## columns, such as in dense HDF5 datasets, and has a fraction of nonzero
## values of at most 'maxparloaddensity' in .memory_model_defaults, estimated
## from a sample of its chunks, see .sample_nonzero_fraction(), the workers
## read it as sparse, see .load_sparse_by_blocks(), and this process converts
## the result into a dense matrix, which, with the parallel reading of the
## chunks, which takes most of the time, and the smaller size of the sparse
## blocks sent back, took less than half the time of reading 100,000 cells of
## single-cell data stored in a dense HDF5 dataset in this process, with six
## workers on macOS. reading data stored by columns, such as in the compressed sparse
## column (CSC) layout of HDF5 files, is instead faster in this process,
## because each column is read without decompressing the values of other
## columns, and so is sending dense blocks back from the workers. within the
## maximum main memory 'maxmem', both the sparse result and the dense matrix
## have to fit, together with the memory held by other objects, see
## .held_mem(), and the one of this process by itself, see
## .block_worker_mem(), and otherwise the data are read in this process
#' @importFrom BiocParallel bpnworkers
#' @importFrom BiocGenerics type
.load_dense <- function(X, maxmem=Inf, BPPARAM=NULL) {
    par <- .parallel_dense_load(X, BPPARAM)
    if (!par)
        return(as.matrix(X))

    density <- attr(par, "density")
    eltbytes <- if (type(X) == "integer") 4 else 8
    sparsesize <- density * as.numeric(nrow(X)) * ncol(X) * (eltbytes + 4)
    if (!is.null(maxmem) && is.finite(maxmem)) {
        densesize <- as.numeric(nrow(X)) * ncol(X) * eltbytes
        pm <- .block_worker_mem(X, 2L, forked=.forked_workers(BPPARAM))
        if (.mem_fraction_R() * maxmem <
            sparsesize + densesize + .held_mem() + pm$workermem)
            return(as.matrix(X))
    }
    res <- .load_sparse_by_blocks(X, maxmem=maxmem, insize=sparsesize,
                                  BPPARAM=BPPARAM, parallelonly=TRUE)
    if (is.null(res))
        return(as.matrix(X))

    as.matrix(res)
}

## whether the on-disk data in 'X' can be loaded into main memory as a dense
## matrix by the parallel workers of 'BPPARAM', see .load_dense(): with at
## least two of them, when 'X' is stored in chunks of several columns and
## has more columns than one chunk, and with a fraction of nonzero values of
## at most 'maxparloaddensity' in .memory_model_defaults, estimated from a
## sample of its chunks, see .sample_nonzero_fraction(), which the attribute
## "density" of the result gives
#' @importFrom BiocParallel bpnworkers
.parallel_dense_load <- function(X, BPPARAM=NULL) {
    if (is.null(BPPARAM) || bpnworkers(BPPARAM) < 2L)
        return(FALSE)
    cd <- .ondisk_chunkdim(X)
    if (is.null(cd) || is.na(cd[2]) || cd[2] < 2 || ncol(X) <= cd[2])
        return(FALSE)

    density <- .sample_nonzero_fraction(X)

    structure(density <= .memory_model("maxparloaddensity"),
              density=density)
}

## with the parallel workers of 'BPPARAM', the sparse data of the GSVA and
## average methods is loaded by them, see .load_sparse_by_blocks(), and the
## data loaded as a dense matrix may be loaded by them, see .load_dense()
.check_sparse_load_input_expr <- function(expr, method, first, last, whdim,
                                          ondisk, verbose, BPPARAM=NULL) {
    ## methods that deal with sparse data, which is loaded in main memory as
    ## sparse, while the other methods load it as a dense matrix
    sparsemethods <- c("GSVA", "average")
    if (!method %in% sparsemethods && is_sparse(expr)) { 
        msg <- paste("Input expression data is sparse, but the {method}",
                     "algorithm does not deal with sparsity in a specific way,",
                     "and data will be converted into a dense matrix format")
        cli_alert_warning(msg)
    }

    if (!is.na(first) || !is.na(last)) {
        if (is.na(first))
            first <- 1
        if (is.na(last))
            last <- if (whdim == 1) nrow(expr) else ncol(expr)
        if (whdim == 1) {
            if (verbose) {
                msg <- "Restricting to rows {first}:{last}"
                cli_alert_info(msg)
            }
            expr <- expr[seq.int(first, last, 1), , drop=FALSE]
        } else if (whdim == 2) {
            if (verbose) {
                msg <- "Restricting to columns {first}:{last}"
                cli_alert_info(msg)
            }
            expr <- expr[, seq.int(first, last, 1), drop=FALSE]
        } else
            cli_abort(c("x"="Invalid internal value for 'whdim' argument."))
    }

    if (is(expr, "DelayedMatrix") && !ondisk) {
        if (verbose)
            cli_alert_info("Loading input expression data into main memory")

        if (method %in% sparsemethods && is_sparse(expr))
            expr <- .load_sparse_by_blocks(expr, maxmem=attr(ondisk, "maxmem"),
                                           insize=attr(ondisk, "insize"),
                                           BPPARAM=BPPARAM)
        else
            expr <- .load_dense(expr, maxmem=attr(ondisk, "maxmem"),
                                BPPARAM=BPPARAM)

        if (ncol(expr) > 10000) ## free up ASAP memory we need not anymore
          invisible(gc())       ## and was allocated when reading a big object
    } 
 
    expr
}

## number of nonzero values of the assay 'assay' of the expression data
## 'exprData', counted when it is in main memory, when its seeds store it, see
## .seeds_nzcount(), or when it has fewer than 2000 columns, and otherwise
## estimated from a sample of its blocks, see .sample_nonzero_fraction().
## with 'nzmin=TRUE', the smallest nonzero value of sparse data is given in
## the attribute "nzmin" of the result, see .nzmin(), found among all its
## nonzero values when they are counted from the data, and otherwise among
## those of a sample of its blocks, as the attribute "nzminexact" tells
#' @importFrom S4Arrays is_sparse
#' @importFrom SparseArray nzcount
#' @importFrom cli cli_alert_info cli_abort
.estimate_nzcount <- function(exprData, assay, verbose, nzmin=FALSE) {
    X <- unwrapData(exprData, assay)
    ## coerce to double to ensure we can deal with numbers larger than 2^31
    nr <- as.numeric(nrow(X))
    nc <- as.numeric(ncol(X))
    nzc <- tot <- nr*nc
    nzm <- NA_real_
    nzmexact <- TRUE
    if (is_sparse(X)) {
        estimated_flag <- FALSE
        if (is(X, "dgCMatrix") || is(X, "SVT_SparseMatrix")) {
            nzc <- as.numeric(nzcount(X))
            if (nzmin)
                nzm <- .nzmin_blocks(X)
        } else if (is(X, "DelayedMatrix")) {
            seedsnzc <- .seeds_nzcount(X)
            if (!is.na(seedsnzc)) {
                nzc <- seedsnzc
                if (nzmin) {
                    nzm <- attr(.sample_nonzero_fraction(X, nzmin=TRUE),
                                "nzmin")
                    nzmexact <- FALSE
                }
            } else if (nc < 2000) {
                Xc <- as(X, "dgCMatrix")
                nzc <- as.numeric(nzcount(Xc))
                if (nzmin)
                    nzm <- .nzmin(Xc)
                rm(Xc)
            } else {
                frac <- .sample_nonzero_fraction(X, nzmin=nzmin)
                nzc <- ceiling(tot * as.vector(frac))
                nzm <- if (nzmin) attr(frac, "nzmin") else NA_real_
                nzmexact <- FALSE
                estimated_flag <- TRUE
            }
        } else
            cli_abort(c("x"="{class(X)} sparse matrix class cannot be handled."))

        if (verbose) {
            estmsg <- ""
            if (estimated_flag)
                estmsg <- " (estimated)"
            cli_alert_info(sprintf("%.0f nonzeros (%s than 2^31) and %.2f%% sparsity%s",
                                   nzc, ifelse(nzc > .Machine$integer.max, "more", "less"),
                                   100 - (100 * nzc / tot), estmsg))
        }
    }
    if (nzmin)
        nzc <- structure(nzc, nzmin=nzm, nzminexact=nzmexact)

    return(nzc)
}

## smallest nonzero value of the matrix or block 'x', ignoring missing values,
## or Inf when it has none
#' @importFrom SparseArray nzvals
.nzmin <- function(x) {
    v <- if (is(x, "dgCMatrix")) x@x
         else if (is(x, "SparseArray")) nzvals(x)
         else x[x != 0]
    v <- v[!is.na(v)]
    if (length(v) == 0L)
        return(Inf)

    min(v)
}

## smallest nonzero value of the sparse matrix 'X' in main memory, see
## .nzmin(), through blocks of columns of the default size of the DelayedArray
## package, so that only the nonzero values of one block are copied at a time
#' @importFrom DelayedArray colAutoGrid read_block
.nzmin_blocks <- function(X) {
    if (is(X, "dgCMatrix"))
        return(.nzmin(X))
    grid <- colAutoGrid(X)

    min(Inf, vapply(seq_along(grid), function(i)
                        .nzmin(read_block(X, grid[[i]], as.sparse=TRUE)),
                    numeric(1)))
}

## number of nonzero values of the DelayedArray object 'x', read from its
## seeds without reading their values, when 'x' only binds seeds that store
## it, such as the HDF5 files of the outputs of gsvaMap() bound by
## gsvaReduce(), possibly setting their dimension names, and NA otherwise,
## e.g., when 'x' subsets its seeds or operates on their values. sampling the
## values of such bound outputs instead opens many of their files for each
## sampled block, which can take minutes on the shared filesystems of clusters
#' @importClassesFrom DelayedArray DelayedArray DelayedSetDimnames DelayedAbind
#' @importClassesFrom HDF5Array H5SparseMatrixSeed
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SparseArray SVT_SparseMatrix
#' @importFrom SparseArray nzcount
.seeds_nzcount <- function(x) {
    if (is(x, "DelayedArray"))
        x <- x@seed
    if (is(x, "DelayedSetDimnames"))
        return(.seeds_nzcount(x@seed))
    if (is(x, "DelayedAbind"))
        return(sum(vapply(x@seeds, .seeds_nzcount, numeric(1))))
    if (is(x, "H5SparseMatrixSeed") || is(x, "dgCMatrix") ||
        is(x, "SVT_SparseMatrix"))
        return(as.numeric(nzcount(x)))

    NA_real_
}

## fraction of nonzero values of the matrix 'X' stored on disk, estimated
## from a sample of up to 'nsample' blocks spread over 'X'. when 'X' is stored
## in chunks, the blocks are chunks, which are read and decompressed at once,
## taken from different rows and columns of chunks, because the fraction of
## nonzero values often differs between rows, e.g., between the genes of
## single-cell data, which may be sorted in the file; otherwise, they are
## blocks of columns with up to 'maxvalues' values each, which data stored by
## columns, such as in CSC layout, read at once, see .h5_csc_create(). with
## 'nzmin=TRUE', the smallest nonzero value of the blocks is given in the
## attribute "nzmin" of the result, see .nzmin()
#' @importFrom S4Arrays ArrayViewport read_block
#' @importFrom SparseArray nzcount
#' @importFrom IRanges IRanges
.sample_nonzero_fraction <- function(X, nsample=32L, maxvalues=1e6,
                                     nzmin=FALSE) {
    dims <- dim(X)
    if (any(dims == 0))
        return(if (nzmin) structure(0, nzmin=Inf) else 0)
    cd <- .ondisk_chunkdim(X)
    if (is.null(cd))
        cd <- c(dims[1L], max(1, floor(maxvalues / as.numeric(dims[1L]))))
    cd <- pmax(1, pmin(cd, dims))
    nchunks <- ceiling(dims / cd)
    k <- seq_len(min(nsample, prod(nchunks)))
    ## the i-th block is in the i-th of the k spaced rows of chunks, and in
    ## the i-th of the k spaced columns of chunks, shifted by half of them
    ## every other block, so that the blocks do not line up
    rowpos <- floor((k - 1) * nchunks[1L] / length(k))
    colpos <- (floor((k - 1) * nchunks[2L] / length(k)) +
               floor(nchunks[2L] / 2) * (k %% 2L)) %% nchunks[2L]
    nz <- tot <- 0
    nzm <- Inf
    for (i in k) {
        vpstart <- c(rowpos[i], colpos[i]) * cd + 1
        vpend <- pmin(vpstart + cd - 1, dims)
        block <- read_block(X, ArrayViewport(dims, IRanges(vpstart, vpend)))
        nz <- nz + as.numeric(nzcount(block))
        tot <- tot + prod(as.numeric(vpend - vpstart + 1))
        if (nzmin)
            nzm <- min(nzm, .nzmin(block))
    }

    if (nzmin) structure(nz / tot, nzmin=nzm) else nz / tot
}
