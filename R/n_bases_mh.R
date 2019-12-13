#' Calculate the number of bases that are homologous to the 3' flanking sequence
#'
#' @description Helper function for extractSigsIndel(). Scans a maximum of 1 indel length from
#' 5' to 3' in the flanking sequence to detect identical bases. Stops prematurely if a
#' non-identical base is found.
#'
#' DSBs can be repaired using 3' microhomology, which can be present relative to either the +
#' or - strand. Therefore, both left and right flanking sequences of the indel need to be
#' considered. When calculating left flanking homology (i.e. homology in the 3' direction for
#' the - strand), the reverse complement of the indel sequence and flanking sequence need to
#' be taken. However, the reverse can be taken for the calculation to save computation.
#'
#' @param indel.seq The indel sequence as a character
#' @param flank.seq The flanking sequence as a character
#'
#' @return An integer stating the number of bases in microhomology
n_bases_mh <- function(indel.seq, flank.seq){
  #indel.sequence = "CTA"
  #flank.sequence = "C"
  
  indel_len <- nchar(indel.seq)
  indel.seq <- unlist(strsplit(indel.seq, ''))
  flank.seq <- unlist(strsplit(flank.seq, ''))[1:indel_len]
  
  n_bases <- 0
  for(i in 1:length(indel.seq)){
    if(indel.seq[i] != flank.seq[i]){
      break
    } else {
      n_bases <- n_bases + 1
    }
  }
  return(n_bases)
}
