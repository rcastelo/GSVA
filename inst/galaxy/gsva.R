options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
})

# Explicitly set locale to C for robust file reading
old_locale <- Sys.getlocale("LC_ALL")
Sys.setlocale("LC_ALL", "C")

suppressPackageStartupMessages({
    library(GSVA)
    library(optparse)
})

option_list <- list(
    make_option(c("-e", "--expression_file"), type = "character", help = "Path to expression data file"),
    make_option(c("-g", "--gene_sets_file"), type = "character", help = "Path to gene sets file"),
    make_option(c("-o", "--output_file"), type = "character", help = "Path to output file"),
    make_option(c("-m", "--method"), type = "character", default = "gsva", help = "GSVA method to use (gsva or ssgsea)"),
    # gsvaParam specific
    make_option(c("-k", "--kcdf"), type = "character", default = "Gaussian", help = "Kernel density estimation function ('Gaussian' or 'Poisson')"),
    make_option(c("-t", "--tau"), type = "numeric", default = 1.0, help = "Exponent for weighting the random walk"),
    make_option(c("-d", "--max_diff"), type = "logical", default = TRUE, help = "Whether to calculate enrichment score as the difference between positive and negative random walk deviations or as the maximum deviation from zero"),
    # ssgseaParam specific
    make_option(c("-a", "--alpha"), type = "numeric", default = 0.25, help = "Exponent for weighting the tail in the random walk for ssGSEA"),
    make_option(c("-n", "--normalize"), type = "logical", default = TRUE, help = "Whether to normalize the scores by the absolute difference between the minimum and maximum for ssGSEA")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

expression_file <- args$expression_file
gene_sets_file <- args$gene_sets_file
output_file <- args$output_file
method <- args$method

# Read expression data with explicit parameters for robustness
expr <- read.table(expression_file, header = TRUE, row.names = 1, sep = ",", check.names = FALSE)

# DEBUG: Print dimensions of expr after reading
print(paste("Dimensions of expr:", paste(dim(expr), collapse = "x")))

# Read gene sets
genesets <- readGMT(gene_sets_file)

# Create parameter object and run GSVA
if (method == "gsva") {
    gsvapar <- gsvaParam(as.matrix(expr), genesets, kcdf = args$kcdf, tau = args$tau, maxDiff = args$max_diff)
    es <- gsva(gsvapar)
} else if (method == "ssgsea") {
    ssgseapar <- ssgseaParam(as.matrix(expr), genesets, alpha = args$alpha, normalize = args$normalize)
    es <- gsva(ssgseapar)
} else {
    stop("Invalid GSVA method specified. Choose 'gsva' or 'ssgsea'.")
}

# Write GSVA results
write.table(es, output_file, sep = "\t", quote = FALSE, col.names = NA)

# Revert locale settings
Sys.setlocale("LC_ALL", old_locale)