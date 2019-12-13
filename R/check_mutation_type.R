#' Check the mutation type argument
#' 
#' A function to check if the mutation type argument contains real 
#' mutation types which can be processed by this package
#'
#' @param type (Optional) Character vector stating the mutation types to be processed.
#' If no argument is given, then the function will return 'snv'
#' @return Vector of all possible mutation types from 'type'
#'
#' @examples 
#' ## Check the combination of mutation types
#' type = c("snv","dbs")
#' check_mutation_type(type)
#'
#' @export

check_mutation_type <- function(type)
{
  # Default type is "snv"
  if (missing(type)) {type = c("snv")}
  
  # Translate type to lower case
  type = unname(sapply(type, function(m) tolower(m)))
  
  if (any(type == "all")) { type = c("snv","dbs","indel") }
  else if (!all(type %in% c("snv","dbs","indel"))) 
  {
    stop("One or more mutation types given are unknown")
  }
  
  return(type)
}
