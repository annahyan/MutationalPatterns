#' Make mutation count matrix of 96 trinucleotides with
#' strand information
#'
#' Make a mutation count matrix with 192 features: 96 trinucleotides and 2 strands,
#' these can be transcription or replication strand
#'
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @param ranges GRanges object with the genomic ranges of:
#' 1. (transcription mode) the gene bodies with strand (+/-) information, or
#' 2. (replication mode) the replication strand with 'strand_info' metadata
#' @param mode "transcription" or "replication", default = "transcription"
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#'
#' @return 192 mutation count matrix (96 X 2 strands)
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
#' ## Transcription strand analysis:
#' ## You can obtain the known genes from the UCSC hg19 dataset using
#' ## Bioconductor:
#' # BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19,
#'   mode = "transcription"
#' )
#'
#' ## You can also use a longer context
#' mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19,
#'   mode = "transcription", extension = 2
#' )
#'
#' ## Replication strand analysis:
#' ## Read example bed file with replication direction annotation
#' repli_file <- system.file("extdata/ReplicationDirectionRegions.bed",
#'   package = "MutationalPatterns"
#' )
#' repli_strand <- read.table(repli_file, header = TRUE)
#' repli_strand_granges <- GRanges(
#'   seqnames = repli_strand$Chr,
#'   ranges = IRanges(
#'     start = repli_strand$Start + 1,
#'     end = repli_strand$Stop
#'   ),
#'   strand_info = repli_strand$Class
#' )
#' ## UCSC seqlevelsstyle
#' seqlevelsStyle(repli_strand_granges) <- "UCSC"
#' # The levels determine the order in which the features
#' # will be countend and plotted in the downstream analyses
#' # You can specify your preferred order of the levels:
#' repli_strand_granges$strand_info <- factor(
#'   repli_strand_granges$strand_info,
#'   levels = c("left", "right")
#' )
#'
#' mut_mat_s_rep <- mut_matrix_stranded(grl, ref_genome, repli_strand_granges,
#'   mode = "replication"
#' )
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_matrix}},
#' \code{\link{mut_strand}}
#'
#' @export

mut_matrix_stranded <- function(vcf_list,
                                ref_genome,
                                ranges,
                                mode = "transcription",
                                extension = 1) {


  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }
  if (inherits(vcf_list, "CompressedGRangesList")) {
    gr_sizes <- S4Vectors::elementNROWS(vcf_list)
    gr <- unlist(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- vcf_list
    gr_sizes <- length(gr)
    names(gr_sizes) <- "My_sample"
  } else {
    .not_gr_or_grl(vcf_list)
  }
  strand <- mut_strand(gr, ranges, mode = mode)
  type_context <- type_context(gr, ref_genome, extension)
  mut_mat <- mut_192_occurrences(type_context, strand, gr_sizes)
  return(mut_mat)
}
