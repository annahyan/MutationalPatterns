#' Make mutation count matrix of 96 trinucleotides
#'
#' @description Make 96 trinucleotide mutation count matrix
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @return 96 mutation count matrix
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
#'
#' ## Construct a mutation matrix with a larger context.
#' ## This is most usefull when you have many mutations per sample.
#' mut_mat_extended <- mut_matrix(vcf_list = grl, ref_genome = ref_genome, extension = 2)
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export
mut_matrix <- function(vcf_list, ref_genome, extension = 1) {

  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }

  # Determine nr mutations per sample
  if (inherits(vcf_list, "CompressedGRangesList")) {
    gr_sizes <- S4Vectors::elementNROWS(vcf_list)
    gr <- BiocGenerics::unlist(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- vcf_list
    gr_sizes <- length(gr)
    names(gr_sizes) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  # Determine type and context of all mutations
  type_context <- type_context(gr, ref_genome, extension)

  # Count the type and context to create the mut_mat
  mut_mat <- mut_96_occurrences(type_context, gr_sizes)
  return(mut_mat)
}
