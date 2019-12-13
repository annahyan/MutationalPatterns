#' Compute strand mutation count vector
#' 
#' Compute strand mutation count vector for all mutation types
#'  
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand, 
#' and "-" for unknown
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' 
#' @noRd
#' @return A list of vectors, where each mutation type has a 
#' vector with mutation counts for two strand

mut_strand_occurrences = function(type_context, strand, type)
{
  # Check mutation type
  type = check_mutation_type(type)
  
  # If type_context is not a list of mutation types, try to find
  # which type is present
  if (all(names(type_context) %in% c("context", "types")))
  {
    context_96 = do.call(rbind, strsplit(TRIPLETS_96,""))[,c(1,3,7)]
    context_96 = paste0(context_96[,1], context_96[,2], context_96[,3])
    if (all(unique(type_context$context) %in% context_96) &
        !is.null(unique(type_context$context)))
    {
      type_context = list("snv"=list("types"=type_context$types,
                                     "context"=type_context$context))
      type = "snv"
    }
    else if (all(unique(type_context$types) %in% DBS) &
             !is.null(unique(type_context$types)))
    {
      type_context = list("dbs"=list("types"=type_context$types,
                                     "context"=type_context$context))
      type = "dbs"
    }
    else if (all(unique(type_context$context) %in% INDEL_CONTEXT) &
             !is.null(unique(type_context$context)))
    {
      type_context = list("indel"=list("types"=type_context$types,
                                       "context"=type_context$context))
      type = "indel"
    }
  }
  
  # get possible strand values
  if (class(strand) != "list")
  {
    values = levels(strand)
    
    idx1 = which(strand == values[1])
    idx2 = which(strand == values[2])
  } else 
  {
    values = levels(strand[[type]])
    
    idx1 = which(strand[[type]] == values[1])
    idx2 = which(strand[[type]] == values[2])
  }
  
  # get type context for both vcf subsets
  if (length(strand) == 0)
  {
    type_context_1 = list("types"=NULL,"context"=NULL)
    type_context_2 = list("types"=NULL,"context"=NULL)
  } else 
  {
    type_context_1 = lapply(type_context[[type]], function(x) x[idx1])
    type_context_2 = lapply(type_context[[type]], function(x) x[idx2])
  }
  
  if (type == "snv")
  {
    # make 96-trinucleotide count vector per set
    vector1 = mut_occurrences(type_context_1, type = type)
    vector2 = mut_occurrences(type_context_2, type = type)
    
    # add names
    names_1 = paste(TRIPLETS_96, values[1], sep = "-")
    names_2 = paste(TRIPLETS_96, values[2], sep = "-")
  } else if (type == "dbs")
  {
    # make DBS count vector per set
    vector1 = mut_occurrences(type_context_1, type = type)
    vector2 = mut_occurrences(type_context_2, type = type)
    
    # add names
    names_1 = paste(DBS, values[1], sep = "-")
    names_2 = paste(DBS, values[2], sep = "-")
  } else if (type == "indel")
  { 
    # make indel count vector per set
    vector1 = mut_occurrences(type_context_1, type = type)
    vector2 = mut_occurrences(type_context_2, type = type)
    
    # add names
    names_1 = paste(INDEL_CONTEXT, values[1], sep = "-")
    names_2 = paste(INDEL_CONTEXT, values[2], sep = "-")
  }
  
  # combine vectors in alternating fashion
  # if one strand is not present, make a vector of zeros with the length
  # of the different mutations
  if (length(vector1) == 0)
  {
    vector1 = rep(0, length(names_1))
    names(vector1) = names_1
  } 
  if (length(vector2) == 0)
  {
    vector2 = rep(0, length(names_2))
    names(vector2) = names_2
  }
  
  vector = c(rbind(vector1, vector2))
  names = c(rbind(names_1, names_2))
  names(vector) = names

  return(vector)
}