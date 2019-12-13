#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving the nonnegative least-squares
#' constraints problem.
#' 
#' @param mut_matrix Named list of count matrices 
#' @param signatures Named list of signature matrices (number of mutational features
#' for each signature matrix must be the same as in the corresponding count matrix)
#' @param cutoff Numeric value of absolute signature contribution. Signatures 
#' greater than or equal to the given value are returned. Default = 0
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrices
#'
#' @importFrom pracma lsqnonneg
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
#'
#' @export

least_squares_error_fitting <- function(mut_matrix, signatures, cutoff = 0)
{
    if (class(mut_matrix) != "list" | isEmpty(names(mut_matrix)) | any(names(mut_matrix) == ""))
      stop("'mut_matrix' is not a list or some elements are not named")
    if (class(signatures) != "list" | isEmpty(names(signatures)) | any(names(signatures) == ""))
      stop("'signatures' is not a list or some elements are not named")
  
    contribution = list()
    reconstructed = list()
    
    for (m in names(mut_matrix))
    {
      if (!(m %in% names(signatures))) 
        stop(sprintf("No signatures found for mutation type %s", m))
      
      # make sure dimensions of input matrix are correct
      if (dim(mut_matrix[[m]])[1] != dim(signatures[[m]])[1])
        stop(paste("Mutation matrix and signatures input have",
                   "different number of mutational features"))
      
      n_features = dim(mut_matrix[[m]])[1]
      n_samples = dim(mut_matrix[[m]])[2]
      n_signatures = dim(signatures[[m]])[2]
      lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
      lsq_reconstructed = matrix(NA, nrow=n_features, ncol=n_samples)
      
      # Process each sample
      for (i in 1:ncol(mut_matrix[[m]]))
      {
        y = mut_matrix[[m]][,i]
        lsq = lsqnonneg(signatures[[m]], y)
        lsq_contribution[,i] = lsq$x
        lsq_reconstructed[,i] = signatures[[m]] %*% as.matrix(lsq$x) 
      }
      
      # Add row and col names
      sample_names = colnames(mut_matrix[[m]])
      signature_names = colnames(signatures[[m]])
      mut_type_names = rownames(signatures[[m]])
      
      colnames(lsq_contribution) = sample_names
      rownames(lsq_contribution) = signature_names
      
      lsq_contribution = lsq_contribution[which(rowSums(lsq_contribution)>=cutoff),
                                          ,drop=FALSE] 
      
      colnames(lsq_reconstructed) = sample_names
      rownames(lsq_reconstructed) = mut_type_names
      
      contribution[[m]] = lsq_contribution
      reconstructed[[m]] = lsq_reconstructed
    }
    
    if (length(contribution) == 1)
    {
      contribution = contribution[[1]]
      reconstructed = reconstructed[[1]]
    }
    
    res = list(contribution, reconstructed)
    names(res) = c("contribution", "reconstructed")
    
    return(res)
}
