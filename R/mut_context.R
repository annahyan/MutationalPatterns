#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome. The function is also able to extract
#' the indel context from the COSMIC database or the predefined database given
#' by this package
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @return Character vector with the context of the base substitutions
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings getSeq
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
#' mut_context <- mut_context(vcfs[[1]], ref_genome, type)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_context = function(vcf, ref_genome, type) 
{
    # Make sure that the chromosome names are compatible with each other.
    if (!(all(seqlevels(vcf) %in% seqlevels(get(ref_genome)))))
        stop(paste( "The chromosome names (seqlevels) of the VCF and the",
                    "reference genome object do not match. Use the",
                    "'seqlevelsStyle()' function to rename chromosome",
                    "names.") )
  
    # Check the given mutation type
    type = check_mutation_type(type)
  
    # Double base substitutions have no context
    if ("dbs" %in% type)
        stop("Extract context is not available for double base substitutions")
  
    contexts = list()
    for (m in type)
    {
      if (m == "snv")
      {
        # Select the SNVs from the vcf
        input_vcf <- vcf[nchar(as.character(vcf$REF))==1 & nchar(as.character(unlist(vcf$ALT)))==1]
        
        # Get the bases on position -1 and +1 from the snv
        vcf_context = as.character(getSeq(get(ref_genome),
                                          seqnames(input_vcf),
                                          start(input_vcf) - 1,
                                          end(input_vcf) + 1))
        contexts[[m]] = vcf_context
      } else if (m == "indel")
      {
        ref_len = nchar(as.character(vcf$REF))
        alt_len = nchar(as.character(unlist(vcf$ALT)))
        
        # Select the indels from the vcf
        input_vcf <- vcf[ref_len != alt_len & (ref_len == 1 | alt_len == 1),]
        
        # Make a bed table to be used in extract_indels()
        bed <- data.frame("chrom"=as.character(seqnames(input_vcf)),
                          "pos"=start(input_vcf),
                          "ref"=as.character(input_vcf$REF),
                          "alt"=as.character(unlist(input_vcf$ALT)))
        
        # Default indel context is "cosmic"
        vcf_context = extract_indels(bed, context.database=INDEL,
                                     ref.genome = ref_genome)
        
        contexts[[m]] = vcf_context
      }
    }
    
    # Return a vector when there is only 1 mutation type
    if (length(names(contexts)) == 1)
      contexts = contexts[[1]]
    
    return(contexts)
}

##
## Deprecated variants
##

#'
#' This function has been removed.  Use 'mut_context' instead.
#'
#' @param vcf        A GRanges object
#' @param ref_genome The reference genome
#'
#' @return Character vector with the context of the base substitutions
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
#' mut_context <- mut_context(vcfs[[1]], ref_genome)
#'
#' @seealso
#' \code{\link{mut_context}}
#'
#' @export

mutation_context <- function(vcf, ref_genome)
{
  .Defunct("mut_context", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_context' instead."))
}
