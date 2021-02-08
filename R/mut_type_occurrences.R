#' Count the occurrences of each base substitution type
#'
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @return data.frame with counts of each base substitution type for
#' each sample in vcf_list
#'
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @importFrom magrittr %>%
#'
#' @export

mut_type_occurrences <- function(vcf_list, ref_genome) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  count <- `C>A` <- `C>G` <- `C>T` <- `T>A` <- NULL
  `T>C` <- `T>G` <- `C>T at CpG` <- `C>T other` <- NULL

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
  type_context <- type_context(gr, ref_genome)
  types <- type_context$types

  # Split C>T in C>T other and C>T at CpG
  C_T_i <- types == "C>T"
  types[C_T_i] <- "C>T other"
  CpG_f <- !is.na(BiocGenerics::match(
    type_context$context[C_T_i],
    c("ACG", "CCG", "TCG", "GCG")
  ))
  types[C_T_i][CpG_f] <- "C>T at CpG"
  types <- factor(types, levels = c(
    "C>A", "C>G", "T>A", "T>C", "T>G",
    "C>T at CpG", "C>T other"
  ))

  # Create vector describing the sample of each variant
  sample_vector <- rep(names(gr_sizes), gr_sizes) %>%
    factor(levels = names(gr_sizes))

  # Count per sample then widen into dataframe.
  df <- tibble::tibble("types" = types, "sample" = sample_vector) %>%
    dplyr::group_by(types, sample, .drop = FALSE) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = types, values_from = count) %>%
    dplyr::mutate(`C>T` = `C>T at CpG` + `C>T other`) %>%
    dplyr::select(
      sample, `C>A`, `C>G`, `C>T`, `T>A`,
      `T>C`, `T>G`, `C>T at CpG`, `C>T other`
    ) %>%
    tibble::column_to_rownames("sample")

  return(df)
}
