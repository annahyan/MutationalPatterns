#' Extract mutational signatures from mutation matrix using NMF
#'
#' Decomposes count matrix of mutation type into signatures and contribution of
#' those signatures to the spectra of the samples/vcf files.
#'
#' @param mut_matrix Named list with mutation matrix of single or double substitution and/or indels.
#' Optional to give a matrix object as mutation count matrix. When two mutation types are given in 
#' the same matrix, the data is treated as "combined" and combined signatures will be returned.
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param rank Number of signatures to extract
#' @param nrun (Optional) Number of iterations, default = 200
#' @param method (Optional) Character stating how to use the data. method = "split" will give 
#' results for each mutation type seperately, whereas method = "combine" will give 
#' combined signatures.\cr
#' Default is "split"
#' @return Named list of mutation matrix, signatures and signature contribution for
#' each mutation type (method = "split") or all mutation types together (method = "combine")
#'
#' @import NMF
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## This function is computational intensive.
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

extract_signatures = function(mut_matrix, type, rank, nrun = 200, method = "split")
{
    # Check if mutation matrix is not empty
    if(all(isEmpty(mut_matrix)))
    {
      stop("Provide a named list for 'mut_matrix' with at least one mutation type")
    }
  
    # Check mutation type
    type = check_mutation_type(type)
    
    # If "mut_matrix" is a matrix, then search for the right mutation type. If not found
    # mutation matrix is treated as combination of mutation types
    if (class(mut_matrix) == "matrix")
    {
      if (all(rownames(mut_matrix) %in% TRIPLETS_96)) 
      { 
        mut_matrix = list("snv"=mut_matrix)
        rank = c("snv"=rank)
      }
      else if (all(rownames(mut_matrix) %in% DBS))
      { 
        mut_matrix = list("dbs"=mut_matrix)
        rank = c("dbs"=rank)
      }
      else if (all(rownames(mut_matrix) %in% INDEL_CONTEXT)) 
      { 
        mut_matrix = list("indel"=mut_matrix) 
        rank = c("indel" = rank)
      }
    }
    
    if (class(mut_matrix) == "matrix") 
    { 
      warning("Mutation type of 'mut_matrix' is unknown. Treated as combined mutation types")
      method = "combine"
    }
    
    if (method == "split")
    {
      mut_matrix = mut_matrix[type]
      
      if (length(names(rank)) == 0 & length(rank) == 1)
      {
        rank = rep(rank, length(type))
        names(rank) = type
      } else {
        rank = rank[type]
      }
      
      if (length(mut_matrix) == 0 | length(rank) == 0)
        stop("One or more mutation types not found in count matrices or ranks")
      
      # Check names of list of mutation matrices
      if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
        stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
      }
      
      # Perform NMF for each mutation type
      nmf_res <- lapply(names(mut_matrix), function(m)
      {
        mut = mut_matrix[[m]]
        rank = rank[m]
        # Add a small pseudocount to avoid features with zero counts.
        mut = as.matrix(mut) + 0.0001
        
        # Make sure the rank_range is valid.
        if (!(rank > 0 & rank == round(rank)))
          stop("Rank should be a positive integer")
        
        if (ncol(mut) < max(rank))
          stop(paste( "The rank should be smaller than the number of",
                      "samples in the input matrix.") )
        
        # Calculate NMF
        res = nmf(mut, rank=rank, method="brunet", nrun=nrun, seed=123456)
        
        # Find signatures and contribution of signatures
        signatures = NMF::basis(res)
        contribution = NMF::coef(res)
        
        # Reconstruct mutation matrix
        reconstructed = signatures %*% contribution
        return(list(signatures = signatures,
                    contribution = contribution,
                    reconstructed = reconstructed))
      })
      names(nmf_res) = names(mut_matrix)
      signatures = list()
      contribution = list()
      reconstructed = list()
      
      # Store results of all mutation types, together in the same lists
      for (m in names(mut_matrix))
      {
        signatures[[m]] = nmf_res[[m]][["signatures"]]
        contribution[[m]] = nmf_res[[m]][["contribution"]]
        reconstructed[[m]] = nmf_res[[m]][["reconstructed"]]
      }
      
      # Return a vector when there is only 1 mutation type
      if (length(names(mut_matrix)) == 1)
      {
        signatures = signatures[[1]]
        contribution = contribution[[1]]
        reconstructed = reconstructed[[1]]
      }
      
      nmf_res = list(signatures = signatures,
                     contribution = contribution,
                     reconstructed = reconstructed)
    } else if (method == "combine")
    {
      # Check if mut_matrix is a list or already combined to a matrix
      if (class(mut_matrix) == "list") 
      { 
        mut_matrix = mut_matrix[type]
        rank = rank[type]
        
        if (length(mut_matrix) == 0 | length(rank) == 0)
          stop("One or more mutation types not found in count matrices or ranks")
        
        # Check names of list of mutation matrices
        if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
          stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
        }
        
        mut_matrix = do.call(rbind, mut_matrix) 
      }
      
      # Add a small pseudocount to avoid features with zero counts.
      mut_matrix = mut_matrix + 0.0001
      
      rank = min(rank)
            
      # Make sure the rank_range is valid.
      if (!(rank > 0 & rank == round(rank)))
        stop("Rank should be a positive integer")
      
      if (ncol(mut_matrix) < max(rank))
        stop(paste( "The rank should be smaller than the number of",
                    "samples in the input matrix.") )
      
      # Calculate NMF
      res = nmf(mut_matrix, rank=rank, method="brunet", nrun=nrun, seed=123456)
      
      # Find signatures and contribution of signatures
      signatures = NMF::basis(res)
      contribution = NMF::coef(res)
      
      # Reconstruct mutation matrix
      reconstructed = signatures %*% contribution
      nmf_res = list(signatures = signatures,
                     contribution = contribution,
                     reconstructed = reconstructed)
    }
    return(nmf_res)
}
