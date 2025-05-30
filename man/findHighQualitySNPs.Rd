% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateMappingBiasVcf.R
\name{findHighQualitySNPs}
\alias{findHighQualitySNPs}
\title{Find High Quality SNPs}
\usage{
findHighQualitySNPs(
  mapping.bias.file,
  max.bias = 0.2,
  min.pon = 2,
  triallelic = FALSE,
  vcf.file = NULL,
  genome
)
}
\arguments{
\item{mapping.bias.file}{Generated by \code{\link{calculateMappingBiasVcf}}.}

\item{max.bias}{Maximum mapping bias}

\item{min.pon}{Minimum number of normal samples, useful to get reliable
mapping bias.}

\item{triallelic}{By default, ignore positions with multiple alt alleles.}

\item{vcf.file}{Optional VCF file (for example dbSNP). Needs to be 
bgzip and tabix processed.}

\item{genome}{See \code{readVcf}}
}
\value{
A \code{GRanges} object with mapping bias passing filters. 
If \code{vcf.file} is provided, it will be the variants in the
corresponding file overlapping with the passed variants.
}
\description{
Function to extract high quality SNPs from the mapping bias database.
Useful for generating fingerprinting panels etc.
}
\examples{

normal.panel.vcf <- system.file("extdata", "normalpanel.vcf.gz",
    package = "PureCN")
bias <- calculateMappingBiasVcf(normal.panel.vcf, genome = "h19")

}
\author{
Markus Riester
}
