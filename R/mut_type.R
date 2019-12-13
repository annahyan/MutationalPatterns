#' Retrieve mutation types from a VCF object
#' 
#' A function to extract the mutations from a vcf and translate to
#' the 6 common base substitution types for SNV and to the 78 strand-agnostic
#' mutation types of the DBS from COSMIC.
#' 
#' @param vcf A CollapsedVCF object
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @return List with character vector for each mutation type
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' mut_type(vcfs[[1]], type)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mut_type = function(vcf, type) 
{
    # Check the mutation type argument
    type = check_mutation_type(type)
    
    # Get the mutations from the vcf
    muts = mutations_from_vcf(vcf, type)
    
    if (class(muts) != "list")
    {
      muts = list(muts)
      names(muts) = type
    }
    
    converted = list()
    
    # For snv and dbs, convert the mutations to reverse complement
    # when needed for the context
    for (n in names(muts))
    {
      if (n == "snv")
      {
        types = unlist(muts[[n]])
        types = gsub('G>T', 'C>A', types)
        types = gsub('G>C', 'C>G', types)
        types = gsub('G>A', 'C>T', types)
        types = gsub('A>T', 'T>A', types)
        types = gsub('A>G', 'T>C', types)
        types = gsub('A>C', 'T>G', types)
        
        converted = c(converted,list("snv"=types))
      }
      else if (n == "dbs")
      {
        types = unlist(muts[[n]])
        types = gsub('GT>TG', 'AC>CA', types)
        types = gsub('GT>CG', 'AC>CG', types)
        types = gsub('GT>AG', 'AC>CT', types)
        types = gsub('GT>TC', 'AC>GA', types)
        types = gsub('GT>CC', 'AC>GG', types)
        types = gsub('GT>AC', 'AC>GT', types)
        types = gsub('GT>TA', 'AC>TA', types)
        types = gsub('GT>CA', 'AC>TG', types)
        types = gsub('GT>AA', 'AC>TT', types)
        
        types = gsub('AT>TG', 'AT>CA', types)
        types = gsub('AT>GG', 'AT>CC', types)
        types = gsub('AT>CG', 'AT>CG', types)
        types = gsub('AT>TC', 'AT>GA', types)
        types = gsub('AT>GC', 'AT>GC', types)
        types = gsub('AT>TA', 'AT>TA', types)
        
        types = gsub('GG>TT', 'CC>AA', types)
        types = gsub('GG>CT', 'CC>AG', types)
        types = gsub('GG>AT', 'CC>AT', types)
        types = gsub('GG>TC', 'CC>GA', types)
        types = gsub('GG>CC', 'CC>GG', types)
        types = gsub('GG>AC', 'CC>GT', types)
        types = gsub('GG>TA', 'CC>TA', types)
        types = gsub('GG>CA', 'CC>TG', types)
        types = gsub('GG>AA', 'CC>TT', types)
        
        types = gsub('CG>AT', 'CG>AT', types)
        types = gsub('CG>GC', 'CG>GC', types)
        types = gsub('CG>AC', 'CG>GT', types)
        types = gsub('CG>TA', 'CG>TA', types)
        types = gsub('CG>GA', 'CG>TC', types)
        types = gsub('CG>AA', 'CG>TT', types)
        
        types = gsub('AG>TT', 'CT>AA', types)
        types = gsub('AG>GT', 'CT>AC', types)
        types = gsub('AG>CT', 'CT>AG', types)
        types = gsub('AG>TC', 'CT>GA', types)
        types = gsub('AG>GC', 'CT>GC', types)
        types = gsub('AG>CC', 'CT>GG', types)
        types = gsub('AG>TA', 'CT>TA', types)
        types = gsub('AG>GA', 'CT>TC', types)
        types = gsub('AG>CA', 'CT>TG', types)
        
        types = gsub('GC>TT', 'GC>AA', types)
        types = gsub('GC>CT', 'GC>AG', types)
        types = gsub('GC>AT', 'GC>AT', types)
        types = gsub('GC>TG', 'GC>CA', types)
        types = gsub('GC>CG', 'GC>CG', types)
        types = gsub('GC>TA', 'GC>TA', types)
        
        types = gsub('TA>AT', 'TA>AT', types)
        types = gsub('TA>CG', 'TA>CG', types)
        types = gsub('TA>AG', 'TA>CT', types)
        types = gsub('TA>GC', 'TA>GC', types)
        types = gsub('TA>CC', 'TA>GG', types)
        types = gsub('TA>AC', 'TA>GT', types)
        
        types = gsub('GA>TT', 'TC>AA', types)
        types = gsub('GA>CT', 'TC>AG', types)
        types = gsub('GA>AT', 'TC>AT', types)
        types = gsub('GA>TG', 'TC>CA', types)
        types = gsub('GA>CG', 'TC>CG', types)
        types = gsub('GA>AG', 'TC>CT', types)
        types = gsub('GA>TC', 'TC>GA', types)
        types = gsub('GA>CC', 'TC>GG', types)
        types = gsub('GA>AC', 'TC>GT', types)
        
        types = gsub('CA>TT', 'TG>AA', types)
        types = gsub('CA>GT', 'TG>AC', types)
        types = gsub('CA>AT', 'TG>AT', types)
        types = gsub('CA>TG', 'TG>CA', types)
        types = gsub('CA>GG', 'TG>CC', types)
        types = gsub('CA>AG', 'TG>CT', types)
        types = gsub('CA>TC', 'TG>GA', types)
        types = gsub('CA>GC', 'TG>GC', types)
        types = gsub('CA>AC', 'TG>GT', types)
        
        types = gsub('AA>TT', 'TT>AA', types)
        types = gsub('AA>GT', 'TT>AC', types)
        types = gsub('AA>CT', 'TT>AG', types)
        types = gsub('AA>TG', 'TT>CA', types)
        types = gsub('AA>GG', 'TT>CC', types)
        types = gsub('AA>CG', 'TT>CG', types)
        types = gsub('AA>TC', 'TT>GA', types)
        types = gsub('AA>GC', 'TT>GC', types)
        types = gsub('AA>CC', 'TT>GG', types)
        converted = c(converted,list("dbs"=types))
      }
      else if (n == "indel")
      { 
        converted = c(converted, list("indel"=muts[[n]]))
      }
    }
    
    # Return a vector when there is only 1 mutation type
    if(length(names(converted)) == 1)
      converted = converted[[1]]
    
    return(converted)
}

##
## Deprecated functions
##

#'
#' This function has been renamed.  Use 'mut_type' instead.
#'
#' @param vcf  A GRanges object.
#'
#' @return Character vector with base substitution types
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' mut_type(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#' \code{\link{mut_type}}
#'
#' @export

mutation_types <- function(vcf)
{
    .Defunct("mut_type", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_type' instead."))
}
