#' Compute all pairwise cosine similarities between mutational profiles/signatures
#' 
#' Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices. 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#' 
#' @param mut_matrix1 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param mut_matrix2 96 mutation count matrix (dimensions: 96 mutations X m samples)
#' @return Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
#'
#' @examples
#' ## You can download mutational signatures from the COSMIC database:
#' # sp_url = http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                    package="MutationalPatterns"))
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
#' ## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
#' cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}
#' 
#' @export

cos_sim_matrix = function(mut_matrix1, mut_matrix2)
{
  n_samples1 = ncol(mut_matrix1)
  n_samples2 = ncol(mut_matrix2)
  res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
  
  for(s in 1:n_samples1)
  {
    signal1 = mut_matrix1[,s]
    cos_sim_vector = c()
    for(i in 1:n_samples2)
    {
      signal2 = mut_matrix2[,i]
      cos_sim_vector[i] = cos_sim(signal1, signal2)
    }
    res_matrix[s,] = cos_sim_vector
  }
  rownames(res_matrix) = colnames(mut_matrix1)
  colnames(res_matrix) = colnames(mut_matrix2)
  
  return(res_matrix)
}


#' This function has been renamed to 'cos_sim_matrix'.
#'
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param signatures 96 mutation count matrix (dimensions: 96 mutations X m samples)
#'
#' @return Matrix with pairwise cosine similarities
#'
#' @examples
#' ## You can download mutational signatures from the COSMIC database:
#' # sp_url = http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                    package="MutationalPatterns"))
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
#' ## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
#' cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{cos_sim_matrix}}
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}

explained_by_signatures = function(mut_matrix, signatures)
{
  .Defunct("cos_sim_matrix", package="MutationalPatterns",
           msg=paste("This function has been renamed to",
                     "'cos_sim_matrix'."))
}
