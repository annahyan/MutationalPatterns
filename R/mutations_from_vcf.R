#' Retrieve base substitutions from vcf
#' 
#' A function to extract base substitutions of each position in vcf
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitutions
#' @import GenomicRanges
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' muts = mutations_from_vcf(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mutations_from_vcf = function(vcf) 
{
    ref = as.character(vcf$REF)
    alt = as.character(unlist(vcf$ALT))

    # Allow both uppercase and lowercase column names.
    if (length(ref) == 0)
        ref = as.character(vcf$ref)

    if (length(alt) == 0)
        alt = as.character(vcf$alt)

    # If these columns are still missing, there's nothing we can do.
    if (length(ref) == 0 || length(alt) == 0)
        warning("Some of your data is missing an ALT and/or a REF column.")

    muts = paste(ref, alt, sep=">")
    return(muts)
}
