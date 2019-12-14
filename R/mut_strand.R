#' Find strand of mutations
#' 
#' @details
#' For transcription mode:
#' Definitions of gene bodies with strand (+/-) information should be defined
#' in a GRanges object.
#' 
#' For the base substitutions that are within gene bodies, it is determined whether
#' the "C" or "T" base is on the same strand as the gene definition. (Since
#' by convention we regard base substitutions as C>X or T>X.)
#'
#' Base substitions on the same strand as the gene definitions are considered
#' "untranscribed", and on the opposite strand of gene bodies as "transcribed",
#' since the gene definitions report the coding or sense strand, which is
#' untranscribed.
#'
#' No strand information "-" is returned for base substitutions outside gene
#' bodies, or base substitutions that overlap with more than one gene body on 
#' the same strand.
#' 
#' For replication mode:
#' Replication directions of genomic ranges should be defined in GRanges object.
#' The GRanges object should have a "strand_info" metadata column, 
#' which contains only two different annotations, e.g. "left" and "right", or 
#' "leading" and "lagging". The genomic ranges cannot overlap, to allow only one 
#' annotation per location.
#' 
#' For each base substitution it is determined on which strand it is located.
#' No strand information "-" is returned for base substitutions in unannotated 
#' genomic regions.
#' 
#' With the package we provide an example dataset, see example code.
#' 
#'
#' @param vcf GRanges containing the VCF object
#' @param ranges GRanges object with the genomic ranges of:
#' 1. (transcription mode) the gene bodies with strand (+/-) information, or 
#' 2. (replication mode) the replication strand with 'strand_info' metadata 
#' @param mode "transcription" or "replication", default = "transcription"
#'
#' @return Character vector with transcriptional strand information with
#' length of vcf: "-" for positions outside gene bodies, "U" for
#' untranscribed/sense/coding strand, "T" for
#' transcribed/anti-sense/non-coding strand.
#' 
#' @importFrom GenomicRanges reduce
#'
#' @examples
#' ## For this example we need our variants from the VCF samples, and
#' ## a known genes dataset.  See the 'read_vcfs_as_granges()' example
#' ## for how to load the VCF samples.
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## For transcription strand:
#' ## You can obtain the known genes from the UCSC hg19 dataset using
#' ## Bioconductor:
#' # source("https://bioconductor.org/biocLite.R")
#' # biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' # library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#'
#' ## For this example, we preloaded the data for you:
#' genes_hg19 <- readRDS(system.file("states/genes_hg19.rds",
#'                         package="MutationalPatterns"))
#'
#' mut_strand(vcfs[[1]], genes_hg19, mode = "transcription")
#' 
#' ## For replication strand:
#' ## Read example bed file with replication direction annotation
#' ## Read replistrand data
#' repli_file = system.file("extdata/ReplicationDirectionRegions.bed", 
#'                           package = "MutationalPatterns")
#' repli_strand = read.table(repli_file, header = TRUE)
#' repli_strand_granges = GRanges(seqnames = repli_strand$Chr, 
#'                                ranges = IRanges(start = repli_strand$Start + 1, 
#'                                end = repli_strand$Stop), 
#'                                strand_info = repli_strand$Class)
#' ## UCSC seqlevelsstyle
#' seqlevelsStyle(repli_strand_granges) = "UCSC"
#' 
#' mut_strand(vcfs[[1]], repli_strand_granges, mode = "transcription")
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_strand = function(vcf, ranges, mode = "transcription")
{
  # Transcription mode
  if(mode == "transcription")
  {
    # Reduce gene object to merge gene definitions that overlap on the same strand
    genes = GenomicRanges::reduce(ranges)
    
    # Check consistency of chromosome names.
    if (!(all(seqlevels(vcf) %in% seqlevels(genes))))
      stop(paste( "Chromosome names (seqlevels) of vcf and genes Granges",
                  "object do not match. Use the seqlevelsStyle() function",
                  "to rename chromosome names.") )
    
    # Determine overlap between vcf positions and genes
    overlap = findOverlaps(vcf, genes)
    overlap = as.data.frame(as.matrix(overlap))
    colnames(overlap) = c('vcf_id', 'gene_body_id')
    
    # Remove mutations that overlap with multiple genes and therefore cannot
    # be determined whether they are on transcribed or untranscribed strand
    # duplicated mutations.
    dup_pos = overlap$vcf_id[duplicated(overlap$vcf_id)]
    
    # Index of duplicated mutations
    dup_idx = which(overlap$vcf_id %in% dup_pos)
    
    # Remove all duplicated (non-unique mapping) mutations.
    if (length(dup_idx) > 0)
      overlap = overlap[-dup_idx,]
    
    # Subset of mutations in genes
    vcf_overlap = vcf[overlap$vcf_id]
    
    # Find reference allele of mutations (and strand of reference genome is
    # reported in vcf file).
    ref = vcf_overlap$REF
    
    # Find the strand of C or T (since we regard base substitutions as
    # C>X or T>X) which mutations have ref allele C or T.
    i = which(ref == "C" | ref == "T")
    
    # Store mutation strand info in vector.
    strand_muts = rep(0, nrow(overlap))
    strand_muts[i] = "+"
    strand_muts[-i] = "-"
    
    # Find strand of gene bodies of overlaps.
    strand_genebodies = as.character(strand(genes)[overlap$gene_body_id])
    
    # Find if mut and gene_bodies are on the same strand.
    same_strand = (strand_muts  == strand_genebodies)
    
    # Subset vcf object for both untranscribed and transcribed
    # gene definition represents the untranscribed/sense/coding strand
    # if mutation is on same strand as gene, than its untranscribed.
    U_index = which(same_strand == TRUE)
    
    # If mutation is on different strand than gene, then its transcribed.
    T_index = which(same_strand == FALSE)
    strand = rep(0, nrow(overlap))
    strand[U_index] = "untranscribed"
    strand[T_index] = "transcribed"
    
    # Make vector with all positions in input vcf for positions that do
    # not overlap with gene bodies, report "-".
    strand2 = rep("-", length(vcf))
    strand2[overlap$vcf_id] = strand
    # make factor 
    strand2 = factor(strand2, levels = c("untranscribed", "transcribed", "-"))
  }
  
  # Replication mode
  if(mode == "replication")
  {
    # Check for presence strand_info metadata
    if(is.null(ranges$strand_info))
    {
      stop("GRanges object with genomic regions does not contain 'strand_info' factor as metadata.")
    }
    # Check that only two different annotations 
    if(length(levels(ranges$strand_info)) != 2)
    {
      stop("GRanges object metadata: 'strand_info' factor should contain exactly two different 
           levels, such as 'left' and 'right'.")
    }
    
    # Determine overlap between vcf positions and genomic regions
    overlap = findOverlaps(vcf, ranges)
    overlap = as.data.frame(as.matrix(overlap))
    colnames(overlap) = c('vcf_id', 'region_id')
    
    # remove mutations that overlap with multiple regions
    dup_pos = overlap$vcf_id[duplicated(overlap$vcf_id)]
    # Index of duplicated mutations
    
    dup_idx = which(overlap$vcf_id %in% dup_pos)
    # Remove all duplicated (non-unique mapping) mutations
    if (length(dup_idx) > 0)
    {
      overlap = overlap[-dup_idx,]
      warning("Some variants overlap with multiple genomic regions in the GRanges object. 
              These variants are assigned '-', as the strand cannot be determined.
              To avoid this, make sure no genomic regions are overlapping in your GRanges 
              object.")
    }
    
    # get strand info of region
    strand = ranges[overlap$region_id]$strand_info
    # Make vector with all positions in input vcf for positions that do
    # not overlap with gene bodies, report "-"
    strand2 = rep("-", length(vcf))
    strand2[overlap$vcf_id] = as.character(strand)
    # make factor, levels defines by levels in ranges object
    levels = c(levels(ranges$strand_info), "-")
    strand2 = factor(strand2, levels = levels)
  }
 
  return(strand2)
}

##
## Renamed function
##

#'
#' This function has been renamed. Use 'mut_strand' instead.
#'
#' @param vcf   A GRanges object
#' @param genes The genes
#'
#' @return Character vector with transcriptional strand information
#'
#' @seealso
#' \code{\link{mut_strand}}
#'
#' @export

strand_from_vcf <- function(vcf, genes)
{
  .Defunct("mut_strand", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_strand' instead."))
}
