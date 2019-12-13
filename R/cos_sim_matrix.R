#' Compute all pairwise cosine similarities between mutational profiles/signatures
#' 
#' Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices. 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#' 
#' @param mut_matrix1 Named list of mutation count matrices (dimensions: n mutations X m samples). \cr
#' It is possible to give a matrix. If both 'mut_matrix1' and 'mut_matrix2' are matrices, 
#' cosine similarity between the two is calculated. If one of them is a list, only the cosine 
#' similarity of the mutation type in the matrix is calculated
#' @param mut_matrix2 Named list of mutation count matrices (dimensions: n mutations X m samples)
#' @return Named list of matrices with pairwise cosine similarities 
#' (dimensions: n mutational profiles X m mutational profiles)
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' For named lists, default is 'snv', else default is mutation type of the matrix
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

cos_sim_matrix = function(mut_matrix1, mut_matrix2, type)
{
  # Type is default if not given
  if (missing(type)) {type_default = TRUE}
  else {type_default = FALSE}
  
  # If "mut_matrix1" and "mut_matrix2" are matrices and no type is given
  # then type is unknown and function stopped
  if (class(mut_matrix1) == "matrix" & class(mut_matrix2) == "matrix")
    if (!type_default)
      stop(paste("Mutation type of matrices is unknown, type can not be chosen.",
                 "Remove type argument or give named lists as mutation matrices"))
  
  # Check "type" argument   
  type = check_mutation_type(type)
  
  # If "mut_matrix1" is matrix and "mut_matrix2" is list, then type
  # of "mut_matrix1" is same as type of matrix in "mut_matrix2" with
  # same rownames
  if (class(mut_matrix1) == "matrix" & class(mut_matrix2) == "list")
  { 
    row_list = list()
    for (m in names(mut_matrix2))
    {
      if (all(rownames(mut_matrix2[[m]]) %in% rownames(mut_matrix1))) 
      { 
        row_list[[m]] = mut_matrix1
        mut_matrix1 = row_list
        break 
      }
    }
  } else if (class(mut_matrix1) == "list" & class(mut_matrix2) == "matrix")
  {
    row_list = list()
    for (m in names(mut_matrix1))
    {
      if (all(rownames(mut_matrix1[[m]]) %in% rownames(mut_matrix2))) 
      { 
        row_list[[m]] = mut_matrix2
        mut_matrix2 = row_list
        break 
      }
    }
  }
  # If both are matrices, then set type to default
  if (class(mut_matrix1) == "matrix" & class(mut_matrix2) == "matrix") 
  {
    mut_matrix1 = list("snv"=mut_matrix1)
    mut_matrix2 = list("snv"=mut_matrix2)
  }
  
  res_matrix_list = list()
  
  # Find corresponding types in both arguments
  names_lists = intersect(names(mut_matrix1), names(mut_matrix2))
  if (isEmpty(names_lists)) { stop("No matching mutation type between 'mut_matrix1' and 'mut_matrix2'")}
  if (length(names_lists) < length(names(mut_matrix1)) |
      length(names_lists) < length(names(mut_matrix2)))
    warning("Only cosine similarity matrices for matching mutation types will be computed")
  
  # For each mutation type, calculate the cosine similarity of the matrices
  for (m in names_lists)
  {
    n_samples1 = ncol(mut_matrix1[[m]])
    n_samples2 = ncol(mut_matrix2[[m]])
    res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
    
    for(s in 1:n_samples1)
    {
      signal1 = mut_matrix1[[m]][,s]
      cos_sim_vector = c()
      for(i in 1:n_samples2)
      {
        signal2 = mut_matrix2[[m]][,i]
        cos_sim_vector[i] = cos_sim(signal1, signal2)
      }
      res_matrix[s,] = cos_sim_vector
    }
    rownames(res_matrix) = colnames(mut_matrix1[[m]])
    colnames(res_matrix) = colnames(mut_matrix2[[m]])
    
    res_matrix_list[[m]] = res_matrix
  }

  # If given type not in results, then stop. Else give
  # the results from the given type
  if (any(!(type %in% names(res_matrix_list))) & !type_default)
    stop(paste("One or more given mutation types not found",
               "in matching types between 'mut_matrix1' and 'mut_matrix2'"))
  else
    res_matrix_list = res_matrix_list[type]

  # If there is 1 mutation type as result, give back a matrix
  if (length(res_matrix_list) == 1)
    res_matrix_list = res_matrix_list[[1]]
  return(res_matrix_list)
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
