## This file contains the constants of the memory model of GSVA, which sizes
## the blocks of its calculations within the maximum main memory 'maxmem' and
## estimates the memory that they require, see memory.R and
## block_processing.R. The constants were estimated from the peak memory of
## the calculations measured on single-cell data, on macOS, on Linux and on
## the nodes of a computer cluster, and include a margin; each one is given
## with how it was measured. They can be estimated again with the scripts in
## the directory 'inst/scripts' of the package, and any of them can be
## replaced, without installing GSVA again, with the option
## 'GSVA.memorymodel', a list with the same structure as
## .memory_model_defaults or part of it, e.g.,
## options(GSVA.memorymodel=list(memfractionR=0.65)), see .memory_model().
## These functions are not exported and are only used internally by GSVA.

.memory_model_defaults <- list(
    ## fraction of the maximum main memory available to the allocations of R,
    ## the rest being taken by the memory that R releases but that the
    ## process keeps, and by allocations outside R, e.g., by the HDF5
    ## library; in the steps of GSVA on single-cell data, the peak of memory
    ## allocated by R was 70% to 80% of the peak of memory resident in the
    ## process. on macOS, the parallel workers keep besides, from one round of
    ## blocks to the next, up to about half the memory of a block that R
    ## releases, as measured on single-cell data stored on disk, which this
    ## rest also covers
    memfractionR=0.7,

    ## memory in bytes that each R process takes outside the allocations of
    ## R, such as the memory that R releases but the process keeps, as
    ## measured in forked workers processing blocks of single-cell data stored
    ## on disk: about 0.3 GB on Linux and about 0.45 GB on macOS, where an R
    ## process reading blocks from an HDF5 file takes about 0.65 GB outside R,
    ## including the buffers of 'diskreadmem'; by operating system, as given
    ## by Sys.info(), with 'default' for the others
    procextramem=c(Darwin=0.45 * 1024^3, default=0.3 * 1024^3),

    ## memory in bytes allocated by R in a new R session loading GSVA and the
    ## packages it depends on, such as the one of a job of gsvaMap() or of a
    ## socket worker, about 0.6 GB, which with the memory outside R gives
    ## about 0.9 GB on Linux and about 1.05 GB on macOS
    newprocmem=0.6 * 1024^3,

    ## memory in bytes that each R process takes outside R while it reads
    ## blocks of input data from disk, such as the buffers of the HDF5 library
    ## to decompress the chunks read: about 0.2 GB, measured on single-cell
    ## data stored in chunks of about 1000 x 1000 values
    diskreadmem=0.2 * 1024^3,

    ## memory in bytes per nonzero value of a block of sparse data read from
    ## disk as sparse with the HDF5Array package while it is read, as
    ## measured on single-cell data, see .load_sparse_by_blocks()
    sparsereadbytes=160,

    ## maximum fraction of nonzero values of the on-disk data loaded into main
    ## memory as a dense matrix by parallel workers reading it as sparse, see
    ## .load_dense(): the workers send each block back to this process, which
    ## took about 3 s per GB of the block on macOS, against about 1.5 s per GB
    ## of the dense matrix reading it in this process, so that reading in
    ## parallel is faster up to about a third of nonzero values, from which
    ## this leaves a margin; it concerns time instead of memory, but it was
    ## measured in the same way
    maxparloaddensity=0.25,

    ## memory taken by binding the blocks of columns of a dgCMatrix object,
    ## relative to the size of the result, as measured, unlike binding
    ## SVT_SparseMatrix objects, see .units_per_block()
    dgcassembly=2,

    ## memory taken by realizing each block of columns of the scaled rows of
    ## PLAGE when they are delayed operations on input stored on disk,
    ## besides the size of its dense form, which in total takes about three
    ## times that size, as measured on single-cell data, see plage()
    plagedelayedextra=2,

    ## memory factors of the steps of GSVA, see .step_mem_factors(): the
    ## memory taken by processing each row or column of the input, relative
    ## to its size ('workfactor'), and by its output, relative to that size
    ## ('outfactor'), for sparse or dense input, and for double or integer
    ## input when the output depends on it; 'assembly' gives the memory taken
    ## by the output in main memory, relative to its size, while it is
    ## produced, see .fixed_mem(); 'sparsedisk' the factors relative to the
    ## sparse size of the blocks of sparse data read from disk, see
    ## .ondisk_mem_factors(); and 'forkcopy' the fraction of the memory of the
    ## main R process that each worker forked from it copies, for input in
    ## main memory and read from disk, see .fork_copied_mem()
    stepfactors=list(
        ## the row normalization with CLR reads its input once through blocks
        ## of columns, see .rowStats(), and gives its output in one piece,
        ## which takes twice its size when the input is a dgCMatrix object
        ## converted into an SVT_SparseMatrix object, see .rownorm_clr(); when
        ## the input is stored on disk, its output is a delayed operation that
        ## takes no memory
        rownormclr=list(workfactor=c(sparse=2, dense=5),
                        outfactor=c(double=1, intsparse=1.5, intdense=2),
                        assembly=c(dgc=2, other=1)),
        ## each worker forked from the main R process copies its whole
        ## memory when calculating the row ECDFs of sparse input held in main
        ## memory other than a dgCMatrix object, about 3 GB in each of 4
        ## workers with 100,000 cells, and none of it with the input read
        ## from disk
        rownorm=list(workfactor=c(sparse=6, dense=4),
                     outfactor=c(double=1, intsparse=1.5, intdense=2),
                     forkcopy=c(inmemory=1, ondisk=0)),
        ## the row statistics of each block of columns, four or five doubles
        ## per row, are small and reduced at the end, see .rowStats(); dense
        ## blocks take a copy of the block, the mask of missing values, the
        ## logarithm of their values and temporary matrices, about five times
        ## their size, as measured
        rowstats=list(workfactor=c(sparse=2, dense=5), outfactor=0),
        ## the ranks of a dgCMatrix object are double, as large as it, and
        ## binding their blocks of columns takes twice their size, see
        ## .cbind_blocks(). blocks of sparse data read from disk, whose size is
        ## the one of their dense form, see .units_per_block(), take instead
        ## about 9 to 11 times their sparse size, reading, ranking and writing
        ## them, and up to 15 times with blocks of about a thousand columns,
        ## as measured on single-cell data, so that narrow blocks fit too.
        ## each worker forked from the main R process copies about 0.17 times
        ## its memory when ranking the columns of 100,000 cells held in main
        ## memory
        colranks=list(workfactor=c(sparse=7, dense=2),
                      outfactor=c(intordgc=1, sparse=0.7, dense=0.5),
                      assembly=c(sparsedgc=2),
                      sparsedisk=list(workfactor=14.3, outfactor=0.7),
                      forkcopy=0.2),
        ## the scores of blocks of sparse ranks read from disk take about two
        ## times their sparse size, besides the scores themselves, as measured
        ## on single-cell data
        scores=list(workfactor=2, outfactor=0,
                    sparsedisk=list(workfactor=2.5, outfactor=0)),
        ## the average method gives its scores directly, while PLAGE, z-score
        ## and ssGSEA give them from an intermediate matrix as large as their
        ## input, of scaled rows or column ranks; the scores of z-score, from
        ## its scaled rows, take as much memory as the ones of the average
        ## method, as the ones of PLAGE, from its scaled rows and the singular
        ## vectors of the gene sets, see plage()
        setscores=list(workfactor=2, outfactor=0),
        ## the products of the rows of each gene set in each block of columns
        ## of the scaled rows of PLAGE, whose memory is given separately, see
        ## plage()
        plagegram=list(workfactor=2, outfactor=0),
        rowscaling=list(workfactor=2, outfactor=1),
        ssgsea=list(workfactor=3, outfactor=1),
        ## the scores of ssGSEA, from its ranks, which take from four to six
        ## times the size of a block of integer ranks: the block read from
        ## disk, its integer copy and the temporary vectors of its columns, as
        ## measured on single-cell data, the larger factors with smaller
        ## blocks
        ssgseascores=list(workfactor=6, outfactor=0)))

## constant 'name' of the memory model, see .memory_model_defaults, as
## replaced by the option 'GSVA.memorymodel', where lists, such as the
## memory factors of the steps, are replaced element by element, see
## modifyList()
#' @importFrom utils modifyList
.memory_model <- function(name) {
    value <- .memory_model_defaults[[name]]
    opt <- getOption("GSVA.memorymodel")[[name]]
    if (is.null(opt))
        return(value)
    if (is.list(value) && is.list(opt))
        return(modifyList(value, opt))

    opt
}

## fraction of the maximum main memory available to the allocations of R, see
## .memory_model_defaults
.mem_fraction_R <- function() {
    .memory_model("memfractionR")
}
