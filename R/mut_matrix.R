#' Make mutation count matrix of 96 trinucleotides 
#'  
#' @description Make 96 trinucleotide mutation count matrix
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param num_cores Number of cores used for parallel processing. If no value
#'                  is given, then the number of available cores is autodetected.
#' @return 96 mutation count matrix
#' @import GenomicRanges
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_matrix = function(vcf_list, ref_genome, num_cores)
{
    df = data.frame()

    if (missing(num_cores))
    {
        num_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
            num_cores <- detectCores()
        else
            num_cores = 1
    }

    rows <- mclapply (as.list(vcf_list), function (vcf)
    {
        type_context = type_context(vcf, ref_genome)
        row = mut_96_occurrences(type_context)
        return(row)
    }, mc.cores = num_cores)

    # Merge the rows into a dataframe.
    for (row in rows)
    {
        if (class (row) == "try-error") stop (row)
        df = rbind (df, row)
    }

    names(df) = names(row)
    row.names(df) = names(vcf_list)

    # transpose
    return(t(df))
}
