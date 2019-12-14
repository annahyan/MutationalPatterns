#' Retrieve context of base substitution types
#' 
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitution types.
#'
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @return Mutation types and context character vectors in a named list
#'
#' @importFrom IRanges reverse
#' 
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' type_context <- type_context(vcfs[[1]], ref_genome)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_context}}
#'
#' @export

type_context = function(vcf, ref_genome)
{
    # Deal with empty GRanges objects.
    if (length (vcf) == 0)
    {
        warning("Detected empty GRanges object.")
        res = list(c(), c())
        names(res) = c("types", "context")
        return(res)
    }

    mut_context = mut_context(vcf, ref_genome)
    muts = mutations_from_vcf(vcf)
    types = mut_type(vcf)

    # find the mutations for which the context needs to be adjusted
    x = which(muts != types)

    # subset mut_context
    y = mut_context[x]

    # Change the context of these mutations to reverse complement
    # of the context
    y = reverse(chartr('ATGC', 'TACG', y))

    # replace subset with reverse complement
    mut_context[x] = y

    # return as named list
    res = list(types, mut_context)
    names(res) = c("types", "context")

    return(res)
}