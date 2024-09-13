#' Adjust tumor vs. normal coverage log ratio for tumor purity and ploidy
#' 
#' This function can be used to adjust the log ratio for tumor purity and
#' ploidy for downstream tools that expect a log2 ratio (for example GISTIC).
#' 
#' 
#' @param ratio Vector of log2 tumor vs normal coverage ratios. 
#' @param purity Purity of sample.
#' @param ploidy Ploidy of sample.
#' @param is.log2 \code{log.ratio} is \code{log2} transformed. 
#' @param min.ratio Minimum (non-log2-transformed) ratio. Set to approx -8
#' \code{log2} adjusted.
#' @return \code{numeric(length(log.ratio))}, \code{log.ratio} adjusted 
#' for \code{purity} and \code{ploidy} 
#' @author Markus Riester
#' @references
#   * Zack et al. (2012), Pan-cancer patterns of somatic copy number alteration 
#'    Nature Biotechnology.
#'  * Toal (2018), https://github.com/lima1/PureCN/issues/40
#' 
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt.gz", 
#'     package = "PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt.gz", 
#'     package = "PureCN")
#' normal <- readCoverageFile(normal.coverage.file)
#' tumor <- readCoverageFile(tumor.coverage.file)
#' log.ratio <- calculateLogRatio(normal, tumor)
#' log.ratio.adjusted <- adjustLogRatio(log.ratio, 0.65, 1.73)
#' 
#' @export adjustLogRatio
adjustLogRatio <- function(ratio, purity, ploidy, is.log2 = TRUE, min.ratio = 2^-8) {
    if (is.log2) ratio <- 2^ratio
    adjusted <- (purity * ploidy * ratio + 2 * (1 - purity) * ratio - 2 * (1 - purity)) / (purity * ploidy)
    adjusted <- pmax(min.ratio, adjusted)
    if (is.log2) adjusted <- log2(adjusted)
    return(adjusted)
}

