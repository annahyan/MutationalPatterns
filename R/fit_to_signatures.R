#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving either the nonnegative least-squares
#' constraints problem or the golden ratio search problem.
#' 
#' @param mut_matrix Named list of count matrices 
#' @param signatures Named list of signature matrices (number of mutational features
#' for each signature matrix must be the same as in the corresponding count matrix)
#' @param cutoff (Optional) Numeric value of absolute signature contribution. Signatures 
#' greater than or equal to the given value are returned. Default = 0
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param method (Optional) Character to select the method used to fit the signatures:
#' \itemize{
#'   \item{'least-squares'} {Solve the nonnegative least-squares constraints problem}
#'   \item{'golden-ratio-search'} {Perform gold ratio search algortihm from deconstructSigs package}
#' } 
#' Default is 'least-squares'
#' @param ... Other arguments passed to whichSignatures() when using method = 'golden-ratio-search'
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## You can download the signatures from the COSMIC website:
#' # http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#' 
#' ## Match the order to MutationalPatterns standard of mutation matrix
#' order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
#' ## Reorder cancer signatures dataframe
#' cancer_signatures = cancer_signatures[order,]
#' ## Use trinucletiode changes names as row.names
#' ## row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
#' ## Keep only 96 contributions of the signatures in matrix
#' cancer_signatures = as.matrix(cancer_signatures[,4:33])
#' ## Rename signatures to number only
#' colnames(cancer_signatures) = as.character(1:30)
#'
#'
#' ## Perform the fitting
#' fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#' \code{\link{least_squares_error_fitting}}
#' \code{\link{golden_ratio_search_fitting}}
#'
#' @export

fit_to_signatures = function(mut_matrix, signatures, type, cutoff, method = "least-squares", ...)
{
    # Check mutation type argument
    if (missing(type)) { type_default = TRUE }
    else { type_default = FALSE }
    type = check_mutation_type(type)
    
    # If signature object is a matrix, then look at "mut_matrix" for mutation type
    if (class(signatures) == "matrix")
    { 
      if (class(mut_matrix) == "matrix") 
      { 
        signatures = list("snv"=signatures) 
        mut_matrix = list("snv"=mut_matrix)
      }
      else 
      {
        signatures_list = list()
        for (m in names(mut_matrix))
        {
          if (all(rownames(mut_matrix[[m]]) %in% rownames(signatures))) 
          { 
            signatures_list[[m]] = signatures
            signatures = signatures_list
            break 
          }
        }
      }
      
    # If count matrix object is a matrix, then look at "signatures" for mutation type
    } else if (class(signatures) == "list")
    {
      if (class(mut_matrix) == "matrix") 
      {
        mut_list = list()
        for (m in names(signatures))
        {
          if (all(rownames(signatures[[m]]) %in% rownames(mut_matrix))) 
          { 
            mut_list[[m]] = mut_matrix
            mut_matrix = mut_list
            break 
          }
        }
        
        type = names(mut_matrix)
      }
    }
    
    if (class(mut_matrix) != "list")
      stop(paste("No list is given for 'mut_matrix' and mutation type",
                 "could not be found in signature list"))
    
    # Get the mutation types asked for
    if (!type_default & !(all(type %in% names(mut_matrix))))
      stop("One or more mutation types are not found in count matrices")
    mut_matrix = mut_matrix[intersect(type, names(mut_matrix))]
    signatures = signatures[intersect(type, names(signatures))]
    
    # Extra arguments for whichSignatures()
    dots = list(...)
    
    # Solve the least squares error problem
    if (method == "least-squares")
    {
      if (missing(cutoff)) { cutoff = 0 }
      res = least_squares_error_fitting(mut_matrix, signatures, cutoff)
      
    # Solve the golden ratio search problem
    } else if (method == "golden-ratio-search")
    {
      # If signature.cutoff is not given, but cutoff is, then use 
      # signature.cutoff = cutoff in whichSignatures()
      if (!("signature.cutoff" %in% names(dots)) & !missing(cutoff)) 
        { res = golden_ratio_search_fitting(mut_matrix, signatures, signature.cutoff = cutoff, ...) }
      else 
        { res = golden_ratio_search_fitting(mut_matrix, signatures, ...) }
    }
    
    return(res)
}
