#' Get DBS context
#'
#' Get the DBS COSMIC context on an GRanges/GRangesList object.
#' It applies the get_dbs_context_gr function to each gr in the input,
#' which works by changing the REF and ALT columns of the GRanges into the COSMIC types.
#'
#' @param vcf_list GRanges/GRangesList
#'
#' @return A version of the GRanges/GRangesList object, with modified REF and ALT columns.
#'
#' @seealso
#' \code{\link{get_mut_type}}, \code{\link{read_vcfs_as_granges}}
#' @family DBS
#'
#' @examples
#' ## Get GRangesList with DBS.
#' ## See 'get_mut_type' or 'read_vcfs_as_granges' for more info on how to do this.
#' dbs_grl <- readRDS(system.file("states/blood_grl_dbs.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Set context DBS
#' get_dbs_context(dbs_grl)
#' @importFrom magrittr %>%
#' @export
#'
get_dbs_context <- function(vcf_list) {
  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }
  
  if (inherits(vcf_list, "CompressedGRangesList")){
    
    # Unlist the GrangesList into a single GRanges object.
    gr <- BiocGenerics::unlist(vcf_list, use.names = FALSE)
    
    # Link samples to mutations
    nr_muts <- S4Vectors::elementNROWS(vcf_list)
    gr$INTERNAL_SAMPLENAME <- rep(names(nr_muts), times = nr_muts)
    
    # Get mutation context
    gr <- .get_dbs_context_gr(gr)
    
    # Split the GRanges back into a list.
    sample_indx <- gr$INTERNAL_SAMPLENAME
    gr$INTERNAL_SAMPLENAME <- NULL
    grl <- S4Vectors::split(gr, sample_indx)
    return(grl)
    
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- .get_dbs_context_gr(vcf_list)
    return(gr)
  } else {
    .not_gr_or_grl(vcf_list)
  }
}


#' Get DBS context on a GRanges
#'
#' Get the DBS COSMIC context on an GRanges object.
#' This is done by changing the REF and ALT columns of the GRanges into the COSMIC types.
#'
#' @param gr GRanges
#'
#' @return A version of the GRanges object, with modified REF and ALT columns.
#'
#' @noRd
#'
#' @importFrom magrittr %>%
#'
.get_dbs_context_gr <- function(gr) {

  # Check that no indels are present.
  .check_no_indels(gr)

  ref <- .get_ref(gr)
  rev_main <- as.vector(ref) %in% c("AA", "GG", "AG", "CA", "GA", "GT")
  ref[rev_main] <- Biostrings::reverseComplement(ref[rev_main])

  alt <- unlist(.get_alt(gr))
  alt[rev_main] <- Biostrings::reverseComplement(alt[rev_main])
  alt_v <- as.vector(alt) # By making this a character vector the default %in% can be used.

  # If the reference is identifcal to its reverse complement then the number of possible alts is reduced to 6.
  # Since A CG>TT is the same as a CG>AA.
  # In these cases the rev complement of the alt has to be taken in some cases.
  rev_alt1 <- which(ref == "TA" & alt_v %in% c("CC", "AC", "AG"))
  rev_alt2 <- which(ref == "AT" & alt_v %in% c("GG", "TG", "TC"))
  rev_alt3 <- which(ref == "GC" & alt_v %in% c("CT", "TT", "TG"))
  rev_alt4 <- which(ref == "CG" & alt_v %in% c("AA", "AC", "GA"))
  rev_alt <- c(rev_alt1, rev_alt2, rev_alt3, rev_alt4)
  alt[rev_alt] <- Biostrings::reverseComplement(alt[rev_alt])
  alt <- alt %>%
    as.vector() %>%
    as.list() %>%
    Biostrings::DNAStringSetList()
  gr$REF <- ref
  gr$ALT <- alt
  return(gr)
}
