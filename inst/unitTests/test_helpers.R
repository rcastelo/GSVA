## GPU test support. GPU blocks are skipped with a visible message when the
## device is unusable; export GSVA_REQUIRE_GPU=1 on GPU-designated runners to
## turn any such skip into a hard failure.

.gsva_gpu_skip_note <- new.env(parent = emptyenv())

.gsva_gpu_reason <- function(state) {
    switch(state,
           unbuilt = "GSVA was built without CUDA support",
           driver = "the CUDA driver is missing or too old",
           nodev = "no CUDA-capable device is detected",
           "the CUDA runtime reported an error")
}

.skip_if_no_gpu <- function() {
    state <- GSVA:::.gsva_cuda_state()
    if (identical(state, "ok"))
        return(FALSE)

    msg <- paste0("Skipping GPU assertions: ", .gsva_gpu_reason(state), ".")
    if (!isTRUE(.gsva_gpu_skip_note$notified)) {
        message(msg)
        .gsva_gpu_skip_note$notified <- TRUE
    }

    if (nzchar(Sys.getenv("GSVA_REQUIRE_GPU")))
        stop(msg, call. = FALSE)

    TRUE
}

.run_with_gpu <- function(expr) {
    if (.skip_if_no_gpu())
        return(invisible(NULL))
    expr
}

.check_gpu_equiv <- function(cpu, gpu, tol = 1e-6) {
    if (inherits(cpu, "ExpressionSet")) {
        cpu <- Biobase::exprs(cpu)
    } else if (inherits(cpu, "SummarizedExperiment")) {
        cpu <- SummarizedExperiment::assay(cpu)
    }

    if (inherits(gpu, "ExpressionSet")) {
        gpu <- Biobase::exprs(gpu)
    } else if (inherits(gpu, "SummarizedExperiment")) {
        gpu <- SummarizedExperiment::assay(gpu)
    }

    checkTrue(identical(dim(cpu), dim(gpu)),
              msg = "GPU and CPU scores differ in dim")
    checkTrue(identical(dimnames(cpu), dimnames(gpu)),
              msg = "GPU and CPU scores differ in dimnames")
    checkTrue(all(is.na(cpu) == is.na(gpu)),
              msg = "GPU and CPU scores differ in which entries are missing")
    checkTrue(all(is.nan(cpu) == is.nan(gpu)),
              msg = "GPU emits NaN where the CPU emits NA, or vice versa")
    checkEqualsNumeric(cpu, gpu, tolerance = tol)
}
