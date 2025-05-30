suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(BiocParallel))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--bam"), action = "store", type = "character",
        default = NULL, help = "Input BAM file"),
    make_option(c("--bai"), action = "store", type = "character",
        default = NULL,
        help = "BAM index file. Only necessary for non-standard file naming."),
    make_option(c("--coverage"), action = "store", type = "character",
        default = NULL,
        help = "Input coverage file (supported file formats are GATK and CNVkit)"),
    make_option(c("--intervals"), action = "store", type = "character", default = NULL,
        help = "Interval file as generated by IntervalFile.R"),
    make_option(c("--keep-duplicates"), action = "store_true",
        default = formals(PureCN::calculateBamCoverageByInterval)$keep.duplicates,
        help = "Count reads marked as duplicates [default %default]"),
    make_option(c("--remove-mapq0"), action = "store_true", default = FALSE,
        help = "Not count reads marked with mapping quality 0 [default %default]"),
    make_option(c("--skip-gc-norm"), action = "store_true", default = FALSE,
        help = "Skips GC-normalization [default %default]"),
    make_option(c("--out-dir"), action = "store", type = "character",
        default = NULL,
        help = "Output directory to which results should be written"),
    make_option(c("--chunks"), action = "store", type = "integer",
        default = formals(PureCN::calculateBamCoverageByInterval)$chunks,
        help = "Split intervals into specified number of chunks to reduce memory usage [default %default]"),
    make_option(c("--parallel"), action = "store_true", default = FALSE,
        help = "Use BiocParallel to calculate coverage in parallel whem --bam is a list of BAM files."),
    make_option(c("--cores"), action = "store", type = "integer", default = 1,
        help = "Number of CPUs to use when --bam is a list of BAM files [default %default]"),
    make_option(c("--seed"), action = "store", type = "integer", 
        default = NULL,
        help = "Seed for random number generator [default %default]"),
    make_option(c("-v", "--version"), action = "store_true", default = FALSE,
        help = "Print PureCN version"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE,
        help = "Overwrite existing files")
)

alias_list <- list(
    "keepduplicates" = "keepduplicates",
    "removemapq0" = "remove-mapq0",
    "skipgcnorm" = "skip-gc-norm",
    "outdir" = "out-dir"
)
replace_alias <- function(x, deprecated = TRUE) {
    idx <- match(x, paste0("--", names(alias_list)))
    if (any(!is.na(idx))) {
        replaced <- paste0("--", alias_list[na.omit(idx)])
        x[!is.na(idx)] <- replaced
        if (deprecated) {
            flog.warn("Deprecated arguments, use %s instead.", paste(replaced, collapse = " "))
        }
    }
    return(x)
}
    
opt <- parse_args(OptionParser(option_list = option_list),
    args = replace_alias(commandArgs(trailingOnly = TRUE)),
    convert_hyphens_to_underscores = TRUE)

if (opt$version) {
    message(as.character(packageVersion("PureCN")))
    q(status = 0)
}

if (!is.null(opt$seed)) {
    set.seed(opt$seed)
}

force <- opt$force

bam.file <- opt$bam
index.file <- opt$bai

gatk.coverage <- opt$coverage
interval.file <- opt$intervals

if (is.null(opt$out_dir)) stop("Need --out-dir")

outdir <- normalizePath(opt$out_dir, mustWork = TRUE)
if (file.access(outdir, 2) < 0) stop("Permission denied to write in --outdir.")

interval.file <- normalizePath(interval.file, mustWork = TRUE)

### Calculate coverage from BAM files -----------------------------------------

.checkFileList <- function(file) {
    files <- read.delim(file, as.is = TRUE, header = FALSE)[, 1]
    numExists <- sum(file.exists(files), na.rm = TRUE)
    if (numExists < length(files)) {
        stop("File not exists in file ", file)
    }
    files
}

checkDataTableVersion <- function() {
    if (compareVersion(package.version("data.table"), "1.12.4") < 0) {
        flog.fatal("data.table package is outdated. >= 1.12.4 required")
        q(status = 1)
    }
    return(TRUE)
}

getCoverageBams <- function(bamFiles, indexFiles, outdir, interval.file,
    force = FALSE, keep.duplicates = FALSE, remove_mapq0 = FALSE) {

    bamFiles <- bamFiles
    indexFiles <- indexFiles
    outdir <- outdir
    interval.file <- interval.file
    force <- force

    .getCoverageBam <- function(bam.file, index.file = NULL, outdir,
        interval.file, force) {
        checkDataTableVersion()
        output.file <- file.path(outdir,  gsub(".bam$", "_coverage.txt.gz",
            basename(bam.file)))
        futile.logger::flog.info("Processing %s...", output.file)
        if (!is.null(index.file)) {
            index.file <- normalizePath(index.file, mustWork = TRUE)
            index.file <- sub(".bai$", "", index.file)
        } else if (file.exists(sub("bam$", "bai", bam.file))) {
            index.file <- sub(".bam$", "", bam.file)
        } else {
            index.file <- bam.file
        }
        if (file.exists(output.file) && !force) {
            futile.logger::flog.info("%s exists. Skipping... (--force will overwrite)", output.file)
        } else {
            PureCN::calculateBamCoverageByInterval(bam.file = bam.file,
                interval.file = interval.file, output.file = output.file,
                index.file = index.file, keep.duplicates = keep.duplicates,
                chunks = opt$chunks,
                mapqFilter = if (remove_mapq0) 1 else NA)
        }
        output.file
    }
    BPPARAM <- NULL
    if (!is.null(opt$cores) && opt$cores > 1) {
        suppressPackageStartupMessages(library(BiocParallel))
        BPPARAM <- MulticoreParam(workers = opt$cores)
        flog.info("Using BiocParallel MulticoreParam backend with %s cores.", opt$cores)
    } else if (opt$parallel) {
        suppressPackageStartupMessages(library(BiocParallel))
        BPPARAM <- bpparam()
        flog.info("Using default BiocParallel backend. You can change the default in your ~/.Rprofile file.")
    }

    if (!is.null(BPPARAM) && length(bamFiles) > 1) {
        coverageFiles <- unlist(
            bplapply(seq_along(bamFiles), 
                function(i) .getCoverageBam(bamFiles[i], indexFiles[i], outdir, interval.file, force),
                BPPARAM = BPPARAM)
        )
    } else {
        coverageFiles <-
            sapply(seq_along(bamFiles),
                function(i) .getCoverageBam(bamFiles[i], indexFiles[i], outdir, interval.file, force))
    }

    coverageFiles
}

coverageFiles <- NULL
indexFiles <- NULL

flog.info("Loading PureCN %s...", Biobase::package.version("PureCN"))

suppressPackageStartupMessages(library(PureCN))
debug <- FALSE
if (Sys.getenv("PURECN_DEBUG") != "") {
    flog.threshold("DEBUG")
    debug <- TRUE
}
    
if (!is.null(bam.file)) {
    bam.file <- normalizePath(bam.file, mustWork = TRUE)
    if (grepl(".list$", bam.file)) {
        bamFiles <- .checkFileList(bam.file)
        if (!is.null(index.file)) {
            if (!grepl(".list$", index.file)) {
                stop("list of BAM files requires list of BAI files.")
            }
            indexFiles <- .checkFileList(index.file)
        }
    } else {
        bamFiles <- bam.file
        indexFiles <- index.file
    }
    if (length(bamFiles) != length(indexFiles) && !is.null(indexFiles)) {
        stop("List of BAM files and BAI files of different length.")
    }

    coverageFiles <- getCoverageBams(bamFiles, indexFiles, outdir,
        interval.file, force, opt$keep_duplicates, opt$remove_mapq0)
}

### GC-normalize coverage -----------------------------------------------------
.gcNormalize <- function(gatk.coverage, interval.file, outdir, force) {
    checkDataTableVersion()
    output.file <- file.path(outdir,  gsub(".txt$|.txt.gz$|_interval_summary",
        "_loess.txt.gz", basename(gatk.coverage)))
    outpng.file <- sub("txt.gz$", "png", output.file)
    output.qc.file <- sub(".txt.gz$", "_qc.txt", output.file)

    if (file.exists(output.file) && !force) {
        flog.info("%s exists. Skipping... (--force will overwrite)", output.file)
    } else {
        png(outpng.file, width = 800, height = 800)
        correctCoverageBias(gatk.coverage, interval.file,
            output.file = output.file, output.qc.file = output.qc.file,
                plot.bias = TRUE)
        invisible(dev.off())
   }
}

if (!opt$skip_gc_norm && (!is.null(gatk.coverage) || !is.null(coverageFiles))) {
    # started not from BAMs?
    if (is.null(coverageFiles)) {
        if (grepl(".list$", gatk.coverage)) {
            coverageFiles <- .checkFileList(gatk.coverage)
        } else {
            coverageFiles <- gatk.coverage
        }
    }
    for (gatk.coverage in coverageFiles)
        .gcNormalize(gatk.coverage, interval.file, outdir, force)
}
