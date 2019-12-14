#' Compute 192 mutation count vector
#' 
#' Compute 192 mutation count vector, 96 trinucleotide changes X 2 strands
#'  
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand, 
#' and "-" for unknown
#' 
#' @noRd
#' @return A vector with 192 mutation occurrences and 96 trinucleotides
#' for two strands

mut_192_occurrences = function(type_context, strand)
{
  # get possible strand values
  values = levels(strand)
  
  idx1 = which(strand == values[1])
  idx2 = which(strand == values[2])
  
  # get type context for both vcf subsets
  type_context_1 = lapply(type_context, function(x) x[idx1])
  type_context_2 = lapply(type_context, function(x) x[idx2])
  
  # make 96-trinucleotide count vector per set
  vector1 = mut_96_occurrences(type_context_1)
  vector2 = mut_96_occurrences(type_context_2)
  
  # add names
  names_1 = paste(TRIPLETS_96, values[1], sep = "-")
  names_2 = paste(TRIPLETS_96, values[2], sep = "-")
  
  # combine vectors in alternating fashion
  vector = c(rbind(vector1, vector2))
  names = c(rbind(names_1, names_2))
  names(vector) = names

  return(vector)
}
