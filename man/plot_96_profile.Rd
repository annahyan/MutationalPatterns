% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_96_profile.R
\name{plot_96_profile}
\alias{plot_96_profile}
\title{Plot 96 trinucleotide profile}
\usage{
plot_96_profile(
  mut_matrix,
  colors = NA,
  same_y = TRUE,
  ratios = TRUE,
  condensed = FALSE,
  sample_labels = NA
)
}
\arguments{
\item{mut_matrix}{96 trinucleotide profile matrix
# @param ymax Y axis maximum value, default = 0.2.}

\item{colors}{Optional 6 value color vector.}

\item{same_y}{A boolean describing whether the same y axis should be used for all samples.}

\item{ratios}{Plot ratios of variants or counts. default = TRUE.}

\item{condensed}{More condensed plotting format. Default = F.}

\item{sample_labels}{Labels to include for sample annotation with e.g. TeX.}
}
\value{
96 trinucleotide profile plot
}
\description{
Plot relative contribution of 96 trinucleotides
}
\examples{
## See the 'mut_matrix()' example for how we obtained the
## mutation matrix information:
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

## Plot the 96-profile of three samples
plot_96_profile(mut_mat[, c(1, 4, 7)])

## Plot a condensed profile
plot_96_profile(mut_mat[, c(1, 4, 7)], condensed = TRUE)

## It's also possible to plot signatures, for example signatures
## generated with NMF
## See 'extract_signatures()' on how we obtained these signatures.
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
  package = "MutationalPatterns"
))

## Optionally, provide signature names
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")

## Generate the plot
plot_96_profile(nmf_res$signatures)
}
\seealso{
\code{\link{mut_matrix}},
\code{\link{plot_profile_heatmap}},
\code{\link{plot_river}}
}
