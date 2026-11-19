## This file contains the internal functions that check the availability of a
## CUDA device for calculating GSVA scores on the GPU, resolve the device used
## for the calculations, and give the limits on the size of the gene sets that
## the GPU supports. They are not exported and are only used internally by
## GSVA.

## Check CUDA device availability:
##   "unbuilt" -- package was compiled without CUDA support
##   "ok"      -- at least one CUDA-capable device is usable
##   "driver"  -- CUDA driver missing or too old for this build
##   "nodev"   -- runtime and driver fine, but no device is visible
##   "init"    -- CUDA cannot be initialized in this R process, which is the
##                case of a process forked from one that has already
##                initialized it, such as a worker of a MulticoreParam
##   "error"   -- any other CUDA runtime failure
.gsva_cuda_state <- function() {
    built <- tryCatch({
        getNativeSymbolInfo("gsva_score_genesets_gpu_R", PACKAGE = "GSVA")
        TRUE
    }, error = function(e) FALSE)
    if (!built)
        return("unbuilt")

    .Call("gsva_cuda_state_R", PACKAGE = "GSVA")
}

#' @importFrom cli cli_abort
.check_gpu_device <- function() {
    cuda_state <- .gsva_cuda_state()
    if (cuda_state == "unbuilt") {
        cli_abort(c(
            "x" = paste("GPU mode requested (device = \"gpu\") but this",
                        "package was built without CUDA support."),
            "i" = paste("To use GPU acceleration, reinstall the package",
                        "with the CUDA toolkit available."),
            "i" = "Alternatively, use device = \"cpu\"."
        ))
    } else if (cuda_state == "init") {
        cli_abort(c(
            "x" = paste("GPU mode requested (device = \"gpu\") but CUDA",
                        "cannot be initialized in this R process."),
            "i" = paste("This happens in a process forked from an R session",
                        "that has already used the GPU, such as a worker of",
                        "a MulticoreParam."),
            "i" = paste("Use device = \"cpu\" in forked processes, or",
                        "parallel workers that are new R processes, such",
                        "as those of a SnowParam.")
        ))
    } else if (cuda_state != "ok") {
        reason <- switch(cuda_state,
                         nodev = "no CUDA-capable device is detected",
                         driver = "the CUDA driver is missing or too old",
                         "the CUDA runtime reported an error")
        cli_abort(c(
            "x" = "GPU mode requested (device = \"gpu\") but {reason}.",
            "i" = paste("Check the GPU and CUDA_VISIBLE_DEVICES, or use",
                        "device = \"cpu\".")
        ))
    }
}

## the device calculating the GSVA scores, from the 'device' argument of
## gsva() or gsvaColScores(): "auto" selects the GPU when CUDA can be used in
## this R process, and the data have no missing values, which only the CPU
## handles, see .gsva_score_genesets(); "gpu" gives an error when CUDA cannot
## be used, before any calculation starts
.resolve_device <- function(device, anyna) {
    if (device == "gpu")
        .check_gpu_device()
    else if (device == "auto")
        device <- if (!anyna && .gsva_cuda_state() == "ok") "gpu" else "cpu"

    device
}

.gsva_cuda_max_gset_size <- function() {
    .Call("gsva_cuda_max_gset_size_R", PACKAGE = "GSVA")
}

.gsva_cuda_walk_safe <- function(G, gset_sizes, tau) {
    .Call("gsva_cuda_walk_safe_R", as.integer(G), as.integer(gset_sizes),
          as.double(tau), PACKAGE = "GSVA")
}
