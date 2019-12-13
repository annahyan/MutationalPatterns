#' Retrieve context of mutations
#' 
#' A function to extract the contexts of mutations for all mutation types
#'
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @return Mutation types and context character vectors in a named list
#'
#' @importFrom IRanges reverse
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
#' type_context <- type_context(vcfs[[1]], ref_genome, type)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_context}}
#'
#' @export

type_context = function(vcf, ref_genome, type)
{
    # Check the mutation type argument
    type = check_mutation_type(type)
    
    # Deal with empty GRanges objects.
    if (length (vcf) == 0)
    {
        warning("Detected empty GRanges object.")
        res = list("types"=c(), "context"=c())
        return(res)
    }
    
    res = list()

    for (m in type)
    {
      muts = mutations_from_vcf(vcf, m)
      
      if (isEmpty(muts))
      {
        res[[m]] = list("types"=NULL, "context"=NULL)
        next
      }
      
      types = mut_type(vcf, m)
      
      # if snvs are extracted, convert the base substitutions to the
      # conventional base substitution types
      if (m == "snv")
      {
        mut_context = mut_context(vcf, ref_genome, m)
        
        # find the mutations for which the context needs to be adjusted
        x = which(muts != types)
        
        # subset mut_context
        y = mut_context[x]
        
        # Change the context of these mutations to reverse complement
        # of the context
        y = reverse(chartr('ATGC', 'TACG', y))
        
        # replace subset with reverse complement
        mut_context[x] = y
        
        res[[m]] = list("types"=types, "context"=mut_context)
      } else if (m == "dbs"){
        res[[m]] = list("types"=types)
      } else if (m == "indel"){
        mut_context = mut_context(vcf, ref_genome, m)
        res[[m]] = list("types"=types, "context"=mut_context)
      }
    }
    
    # Return a vector when there is only 1 mutation type
    if (length(res) == 1)
      res = list("types"=res[[1]][["types"]], "context"=res[[1]][["context"]])
    
    return(res)
}