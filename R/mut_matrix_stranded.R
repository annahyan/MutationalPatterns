#' Make mutation count matrix of 96 trinucleotides with 
#' strand information
#'
#' Make a mutation count matrix with 192 features: 96 trinucleotides and 2 strands,
#' these can be transcription or replication strand
#'
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param ranges GRanges object with the genomic ranges of:
#' 1. (transcription mode) the gene bodies with strand (+/-) information, or 
#' 2. (replication mode) the replication strand with 'strand_info' metadata 
#' @param mode (Optional) "transcription" or "replication", default = "transcription"
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param num_cores Number of cores used for parallel processing. If no value
#'                  is given, then the number of available cores is autodetected.
#'
#' @return List of mutation count matrices for each mutation type. Number of different mutations
#' for all types all doubled, because there are 2 strands
#'
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
#' ## Transcription strand analysis:
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
#' mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_hg19, 
#'                                 mode = "transcription")
#' 
#' ## Replication strand analysis:
#' ## Read example bed file with replication direction annotation
#' repli_file = system.file("extdata/ReplicationDirectionRegions.bed", 
#'                           package = "MutationalPatterns")
#' repli_strand = read.table(repli_file, header = TRUE)
#' repli_strand_granges = GRanges(seqnames = repli_strand$Chr, 
#'                                ranges = IRanges(start = repli_strand$Start + 1, 
#'                                end = repli_strand$Stop), 
#'                                strand_info = repli_strand$Class)
#' ## UCSC seqlevelsstyle
#' seqlevelsStyle(repli_strand_granges) = "UCSC"
#' # The levels determine the order in which the features 
#' # will be countend and plotted in the downstream analyses
#' # You can specify your preferred order of the levels:
#' repli_strand_granges$strand_info = factor(repli_strand_granges$strand_info, levels = c("left", "right"))
#' 
#' mut_mat_s_rep = mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
#'                                 mode = "replication")
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_matrix}},
#' \code{\link{mut_strand}}
#'
#' @export

mut_matrix_stranded = function(vcf_list, ref_genome, ranges, mode = "transcription", type, num_cores)
{
  # Check the mutation type
  type = check_mutation_type(type)
  
  df = list()
  
  # Detect number of cores available for parallelization
  if (missing(num_cores))
  {
      num_cores = detectCores()
      if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
          num_cores <- detectCores()
      else
          num_cores = 1
  }

  for (m in type)
  {
    # Transcription mode
    if(mode == "transcription")
    {
      # For each vcf in vcf_list count the mutational features with strand info
      rows <- mclapply (as.list(vcf_list), function (vcf)
      {
        type_context = type_context(vcf, ref_genome, type = m)
        strand = mut_strand(vcf, ranges, type = m, mode = "transcription")
        row = mut_strand_occurrences(type_context, strand, type = m)
        return(row)
      }, mc.cores = num_cores)
      
      df[[m]] = data.frame()
      
      # Combine the rows in one dataframe
      for (row in rows)
        df[[m]] = rbind (df[[m]], row)
    }
  
    # Replication mode
    if(mode == "replication")
    {
      # For each vcf in vcf_list count the 192 features
      rows <- mclapply (as.list(vcf_list), function (vcf)
      {
        type_context = type_context(vcf, ref_genome, type = m)
        if(is.null(type_context$types))
          strand = factor(c(),levels = c("left","right","-"))
        else
          strand = mut_strand(vcf, ranges, type = m, mode = "replication")
        row = mut_strand_occurrences(type_context, strand, type = m)
        return(row)
      }, mc.cores = num_cores, mc.silent = FALSE)
      
      if (isEmpty(do.call(rbind, rows))) { next }
      
      df[[m]] = data.frame()
      
      # Combine the rows in one dataframe
      for (row in rows)
      {
        if (class (row) == "try-error") stop (row)
        df[[m]] = rbind (df[[m]], row)
      }
    }
  
    # Add row names to data.frame
    names(df[[m]]) = names(row)
    row.names(df[[m]]) = names(vcf_list)
    df[[m]] = t(df[[m]])
  }
  
  if (length(df) == 1)
    return(df[[1]])
  else 
    return(df)
}
