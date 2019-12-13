#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by using the golden ratio search
#' algorithm implemented in the \link{deconstructSigs} package (Rosenthal et al. 2016).
#' 
#' @param mut_matrix Named list of count matrices 
#' @param signatures Named list of signature matrices (number of mutational features
#' for each signature matrix must be the same as in the corresponding count matrix)
#' @param ... Other arguments passed to whichSignatures()
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrices
#'
#' @importFrom deconstructSigs whichSignatures
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
#' \code{\link{fit_to_signatures}}
#' \code{\link{whichSignatures}}
#'
#' @export
 
golden_ratio_search_fitting <- function(mut_matrix, signatures, ...)
{
    # Check if mut_matrix and signatures are both lists
    if (class(mut_matrix) != "list" | isEmpty(names(mut_matrix)) | any(names(mut_matrix) == ""))
      stop("'mut_matrix' is not a list or some elements are not named")
    if (class(signatures) != "list" | isEmpty(names(signatures)) | any(names(signatures) == ""))
      stop("'signatures' is not a list or some elements are not named")
  
    mut_matrix_transposed = list()
    signatures_transposed = list()
    contribution = list()
    reconstructed = list()
    unknown = list()
    
    # For each mutation type in mut_matrix
    for (m in names(mut_matrix))
    {
      # Transpose the mutation matrix and signature matrix
      mut_matrix_transposed[[m]] = as.data.frame(t(mut_matrix[[m]]))
      signatures_transposed[[m]] = as.data.frame(t(signatures[[m]]))
      
      # For each sample get results from the golden ratio search
      result = sapply(rownames(mut_matrix_transposed[[m]]), function(n)
        whichSignatures(mut_matrix_transposed[[m]], sample.id = n, 
                        signatures_transposed[[m]], 
                        contexts.needed = TRUE, ...))
      
      # Write results of contribution, reconstructed and unknown as lists 
      # of mutation types
      contribution[[m]] = as.matrix(do.call(cbind, 
                                            lapply(1:ncol(result), function(i)
                                              data.frame(unlist(result[1,i])))))
      colnames(contribution[[m]]) = colnames(result)
      
      reconstructed[[m]] = as.matrix(do.call(cbind, 
                                             lapply(1:ncol(result), function(i) 
                                               data.frame(unlist(result[3,i])))))
      colnames(reconstructed[[m]]) = colnames(result)
      rownames(reconstructed[[m]]) = colnames(result[3][[1]])
      
      unknown[[m]] = as.matrix(do.call(cbind, lapply(1:ncol(result), function(i)
        data.frame(unlist(result[5,i])))))
      colnames(unknown[[m]]) = colnames(result)
    }
    
    if (length(contribution) == 1)
    {
      contribution = contribution[[1]]
      reconstructed = reconstructed[[1]]
      unknown = unknown[[1]]
    }
    
    res = list("contribution"=contribution, 
               "reconstructed"=reconstructed, 
               "unknown"=unknown)
    
    return(res)
}
