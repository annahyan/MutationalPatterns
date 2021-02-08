#' Rename NMF signatures based on previously defined signatures
#'
#' This function renames signatures identified with NMF based on previously defined signatures.
#' If a NMF signature has a cosine similarity with a previously defined signature,
#' that is higher than the cutoff, then this NMF signature will get the name
#' of the previously defined signature. If not the NMF signature will receive a letter based name.
#' For example: SBSA.
#' This only changes the names of signatures, not their actual values.
#' This function can be help with identifying whether signatures found with NMF are already known,
#' which can be useful for interpretation.
#' An extracted signature that is not similar to any previously defined signatures,
#' is not proof of a "novel" signature. The extracted signature might be a
#' combination of known signatures, that could not be split by NMF. This can happen
#' when, for example, too few samples were used for the NMF.
#' 
#' @param nmf_res Named list of mutation matrix, signatures and signature contribution
#' @param signatures A signature matrix
#' @param cutoff Cutoff at which signatures are considered similar. Default: 0.85
#' @param base_name The base part of a letter based signature name. Default: "SBS"
#' @param suffix String. The suffix added to the name of a renamed signature. Default: "-like"
#' 
#' @return A nmf_res with changed signature names
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' rename_nmf_signatures(nmf_res, signatures)
#' 
#' ## You can change or remove the suffix of the renamed signatures.
#' rename_nmf_signatures(nmf_res, signatures, suffix = "")
#'
#' ## You can change how similar the signatures have to be, before they are considered similar.
#' rename_nmf_signatures(nmf_res, signatures, cutoff = 0.95)
#'
#' ## You can also change the base_name of the signatures that end up with a letter name.
#' rename_nmf_signatures(nmf_res, signatures, cutoff = 0.95, base_name = "Signature_")
rename_nmf_signatures <- function(nmf_res, 
                                  signatures, 
                                  cutoff = 0.85, 
                                  base_name = "SBS", 
                                  suffix = "-like") {
  rownames(nmf_res$contribution) <- seq_len(nrow(nmf_res$contribution))
  colnames(nmf_res$signatures) <- seq_len(ncol(nmf_res$signatures))
  sim_matrix <- cos_sim_matrix(signatures, nmf_res$signatures)

  # The number of signatures with no similarity based match.
  j <- 0

  # Iterate over signatures in nmf_res
  for (i in seq_len(ncol(sim_matrix))) {
    sig_sim <- sim_matrix[, i]

    # Check if there is a highly similar signature
    cossim <- max(sig_sim)
    if (cossim > cutoff) {

      # Determine most similar signature
      row <- which.max(sig_sim)
      sig <- names(sig_sim)[row]
      
      # Add suffix
      sig <- paste0(sig, suffix)

      # Change name of nmf_res signature
      rownames(nmf_res$contribution)[i] <- sig
      colnames(nmf_res$signatures)[i] <- sig
    } else {
      # If there is no similar signature, use letters for the signature name.
      j <- j + 1
      rownames(nmf_res$contribution)[i] <- paste0(base_name, LETTERS[j])
      colnames(nmf_res$signatures)[i] <- paste0(base_name, LETTERS[j])
    }
  }

  # Check if there are any signatures that have been given the same name.
  dupli_sigs <- colnames(nmf_res$signatures) %>%
    duplicated() %>%
    any()
  if (dupli_sigs) {
    stop("You have multiple NMF signatures that are linked to the same existing signature.\n
                Please use a lower rank in the NMF or increase the cutoff at which a NMF and \n
                existing signature are considered identical.", call. = FALSE)
  }

  # Return the nmf_res with the updated signature names.
  return(nmf_res)
}
