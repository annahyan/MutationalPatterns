#' Make mutation count matrices of single and double substitutions and indels
#'  
#' @description Make mutation count matrices for 96 trinucleotide single base 
#' substitutions, 78 double base substitutions and indels. Number of indels 
#' depends on the indel context given by the user
#' 
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param indel (Optional) List of mutation matrix and vectors for context, class and color of indels.
#' Color vector must have the same length as the class vector, 
#' because in plotting profiles, each class is represented by one color.\cr\cr
#' It is also possible to give a character for predefined variables
#' and counting mutations with predefined context:
#' \itemize{
#'   \item{"predefined"} {Represents a indel context of 3 classes per 
#'   deletion and indel: "repetitive region", "microhomology" and 
#'   "none" of these two. Indels have lengths 1 to 5+}
#'   \item{"cosmic"} {Represents the indel context according to the
#'   COSMIC database}
#' }
#' Default is "cosmic"
#' @param method (Optional) Character stating how to use the data. 
#' \itemize{
#'   \item{"split":} { Each mutation type has seperate count matrix}
#'   \item{"combine":} { Combined count matrix of all mutation types}
#' }   
#' Default is "split"
#' @param num_cores (Optional) Number of cores used for parallel processing. If no value
#'                  is given, then the number of available cores is autodetected.
#' @param ... Arguments passed to type_context()                  
#'  
#' @return List with mutation count matrices of snv (96 mutations), dbs (78 
#' mutations) and indels (number of mutations depends on chosen context for indels)
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
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome, type)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{type_context}}
#'
#' @export

mut_matrix = function(vcf_list, ref_genome, type, method = "split", num_cores)
{
    
  
    # Check value of method
    if (!(method %in% c("split", "combine"))){ stop("Provide the right value of 'method'. Options are 'split' or 'combine'")}
  
    df = list("snv"=data.frame(), "dbs"=data.frame(), "indel"=data.frame())

    # Set number of cores used for counting mutations
    if (missing(num_cores))
    {
        num_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
            num_cores <- detectCores()
        else
            num_cores = 1
    }
    
    # Check mutation type 
    type = check_mutation_type(type)

    rows <- mclapply (as.list(vcf_list), function (vcf)
    {
        row = list()
        for (m in type)
        {
          # For every present mutation type count the occurrences of mutations
          if (m == "snv") { 
            row[[m]] = mut_occurrences(type_context(vcf, ref_genome, m), type = m) }
          else if (m == "dbs") { row[[m]] = mut_occurrences(type_context(vcf, ref_genome, m), type = m) }
          else if (m == "indel")
          {
            if(INDEL == "custom") 
            {
              # If a custom mutation matrix is given for the indels, select those
              # columns which correspond to samples in the vcfs
              warning("Custom classification of indels used")
              column = which(colnames(indel$matrix) == names(vcf))
              row[[m]] = indel$matrix[,column]
            } else { 
              row[[m]] = mut_occurrences(type_context(vcf, ref_genome, m), type = m) }
          }
        }
        return(row)
    }, mc.cores = num_cores)

    # Merge the rows into a dataframe of each mutation type.
    for (m in names(rows[[1]]))
    {
      rnames = c()
      for (i in 1:length(rows))
      {
        if (class (row) == "try-error") stop (row)
        if (length(rows[[i]][[m]]) > 0)
        {
          rnames = c(rnames, names(vcf_list)[i])
          df[[m]] = rbind(df[[m]], rows[[i]][[m]])
          names(df[[m]]) = names(rows[[i]][[m]])
        }
      }
      
      if (!isEmpty(rnames)) rownames(df[[m]]) = rnames
      df[[m]] <- t(df[[m]])
    }
    
    # If there is one mutation type, return a matrix,
    # else return a list of mutation types
    if (method == "split") 
    { 
      if (length(which(!isEmpty(df))) == 1) 
      {
        return(df[[which(!isEmpty(df))]])
      } else
      {
        return(df[which(!isEmpty(df))])
      }
    } else if (method == "combine") { return(do.call(rbind, df[which(!isEmpty(df))]))}
}
