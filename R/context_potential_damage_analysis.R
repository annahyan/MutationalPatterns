#' Potential damage analysis for the supplied mutational contexts
#'
#' The ratio of possible 'stop gain', 'mismatches', 'synonymous mutations' and
#' 'splice site mutations' is counted per mutational context. This is done for
#' the supplied ENTREZ gene ids. This way it can be determined how damaging a
#' mutational context could be. N gives the total number of possible mutations
#' per context.
#'
#' The function works by first selecting the longest transcript per gene. The
#' coding sequence (cds) of this transcript is then assembled. Next, the
#' function loops over the reference contexts. For each context (and it's
#' reverse complement), all possible mutation locations are determined. Splice
#' site mutations are removed at this stage. It's also determined whether these
#' locations are the first, second or third base of the cds codon (mut loc).
#' Each unique combination of codon and mut loc is then counted. For each
#' combination the reference amino acid and the possible alternative amino acids
#' are determined. By comparing the reference and alternative amino acids, the
#' number of 'stop_gains', 'mismatches' and 'synonymous mutations' is
#' determined. This is then normalized per mutation context.
#' For example, mutations with the ACA context could be located in the third
#' position of a codon like TAC. This might happen 200 times in the supplied
#' genes. This TAC codon could then be mutated in either a TAA, TAG or a TAT.
#' The first two of these options would induce a stop codon, while the third one
#' would be synonymous. By summing up all codons the number of stop_gains',
#' 'mismatches' and 'synonymous mutations' is determined per mutation context.
#' 
#' For mismatches the blosum62 score is also calculated. This is a score based
#' on the BLOSUM62 matrix, that describes how similar two amino acids are. This
#' score is normalized over the total amount of possible mismatches. A lower
#' score means that the amino acids in the mismatches are more dissimilar. More
#' dissimilar amino acids are more likely to have a detrimental effect. 
#' 
#' To identify splice sites, sequences around the splice locations are used
#' instead of the cds. The 2 bases 5' and 2 bases 3' of a splice site are
#' considered to be splice site mutation locations.
#' 
#' @param contexts Vector of mutational contexts to use for the analysis.
#' @param txdb Transcription annotation database
#' @param ref_genome BSgenome reference genome object
#' @param gene_ids Entrez gene ids
#' @param verbose Boolean. Determines whether progress is printed. (Default: FALSE)
#'
#' @return A tibble with the ratio of 'stop gain', 'mismatch', 'synonymous' and 
#' 'splice site' mutations per mutation context.
#' @export
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the
#' ## mutation matrix information:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' contexts <- rownames(mut_mat)
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Load the transcription annotation database
#' ## You can obtain the database from the UCSC hg19 dataset using
#' ## Bioconductor:
#' # BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Here we will use the Entrez Gene IDs from several cancer
#' ## genes. In practice you might want to use larger gene lists,
#' ## but here we only use a few to keep the runtime low.
#' ## In this example we are using:
#' ## TP53, KRAS, NRAS, BRAF, BRCA2
#' gene_ids <- c(7157, 3845, 4893, 673, 675)
#'
#' ## Run the function
#' context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)
#'
#' ## The function can provide updates about its progress.
#' ## This can be usefull when it's running slowly,
#' ## which can happen when you are using many gene_ids.
#' ## To reduce the example runtime, we don't re-run the analysis, but only show the command
#' ## context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids, verbose = TRUE)
#' 
context_potential_damage_analysis <- function(contexts, txdb, ref_genome, gene_ids, verbose = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  context <- n <- ratio <- NULL
  
  # Check dependencies are installed.
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(paste0(
      "Package 'GenomicFeatures' is needed for context_potential_damage_analysis to work. ",
      "Please install it if you want to use this function."
    ), call. = FALSE)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(paste0(
      "Package 'AnnotationDbi' is needed for context_potential_damage_analysis to work. ",
      "Please install it if you want to use this function."
    ), call. = FALSE)
  }

  # Get reference genome
  ref_genome <- BSgenome::getBSgenome(ref_genome)
  
  # Get DNA of exons. This is strand specific, so I don't need to worry about this.
  cds_tx <- .get_cds_ranges(txdb, gene_ids)
  seqs <- .get_cds_sequences(cds_tx, ref_genome)
  
  splice_genelocs <- .get_splice_genelocs(cds_tx)
  
  if (verbose) {
    message("Finished getting the coding sequences.")
  }

  # Get substitution and contexts
  substitution <- stringr::str_replace(contexts, "\\w.*\\[(.*)\\]\\w.*", "\\1")
  l_context <- stringr::str_remove(contexts, "\\[.*")
  r_context <- stringr::str_remove(contexts, ".*\\]")
  ori_bases <- stringr::str_replace(contexts, "\\[(.*)>.*\\]", "\\1")
  ref_base <- stringr::str_remove(substitution, ">.*")
  alt_base <- stringr::str_remove(substitution, ".*>")

  # Group by ori_bases, because the DNA needs to be searched based on them.
  contexts_tb <- tibble::tibble(
    "ori_bases" = ori_bases,
    "ref_base" = ref_base,
    "alt_base" = alt_base,
    "l_context" = l_context,
    "r_context" = r_context
  ) %>%
    dplyr::group_by(ori_bases) %>%
    dplyr::summarise(
      ref_base = ref_base[[1]],
      alt_bases = list(alt_base),
      l_context = l_context[[1]],
      r_context = r_context[[1]],
      .groups = "drop_last"
    ) %>% 
    dplyr::mutate(mut_pos = stringr::str_length(l_context) + 1,
                  rev_mut_pos = stringr::str_length(r_context) + 1)

  # Read in PAM matrix
  blosum62 <- readRDS(system.file(file.path("states", "blosum62.rds"), 
                               package = "MutationalPatterns"))

  # Perform damage analysis per context.
  mismatches <- purrr::map(seq_len(nrow(contexts_tb)), 
                           .single_context_damage_analysis, 
                           contexts_tb, 
                           seqs, 
                           verbose,
                           blosum62,
                           splice_genelocs) %>%
    dplyr::bind_rows() %>% 
    dplyr::mutate(context = factor(context, levels = unique(context)))

  # Perform damage analysis for splice sites
  splice_muts_tb <- .potential_splice_site_damage(contexts_tb, cds_tx, ref_genome)

  #Combine mismatches tibble with splice site tibble
  mismatches <- rbind(mismatches, splice_muts_tb) %>%
    dplyr::arrange(context) %>% 
    dplyr::group_by(context) %>%
    dplyr::mutate(ratio = n / sum(n)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(type, context, n, ratio, blosum62)
  
  return(mismatches)
}


#' Get cds sequences for supplied genes
#'
#' Per gene the longest transcript is used.
#'
#' @param cds_tx GRangesList object containing the cds ranges of genes
#' @param ref_genome BSgenome reference genome object
#'
#' @return DNAStringSet containing the cds sequences
#' @noRd
#'
.get_cds_sequences <- function(cds_tx, ref_genome) {
  
  # Get sequences (per cds per transcript.)
  seqs <- Biostrings::getSeq(ref_genome, cds_tx)

  # Merge cds sequences per transcript
  seqs <- purrr::map(as.list(seqs), function(seq) do.call(c, as.list(seq))) %>%
    Biostrings::DNAStringSet()

  return(seqs)
}

#' Get cds ranges for supplied genes
#'
#' Per gene the longest transcript is used.
#'
#' @param txdb Transcription annotation database
#' @param gene_ids Entrez gene ids
#'
#' @return GRangesList containing the cds ranges of genes
#' @noRd
#'
.get_cds_ranges = function(txdb, gene_ids){
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  GENEID <- tx_size <- TXNAME <- NULL
  
  # Get cds per transcript
  cds_tx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  
  # Get sizes of transcripts
  tx_sizes <- cds_tx %>%
    BiocGenerics::width() %>%
    sum() %>%
    tibble::enframe(name = "TXNAME", value = "tx_size")
  
  # Get gene names belonging to transcripts
  withCallingHandlers(
    { # Supress the returned 1:many mapping message
      gene2txname <- AnnotationDbi::select(txdb, AnnotationDbi::keys(txdb, "GENEID"), 
                                           columns = c("GENEID", "TXNAME"), 
                                           keytype = "GENEID")
    },
    message = function(m) {
      if (grepl(" returned 1:many mapping between keys and columns", conditionMessage(m))) {
        invokeRestart("muffleMessage")
      }
    }
  )
  gene2txname <- gene2txname %>%
    dplyr::filter(GENEID %in% gene_ids) %>%
    dplyr::inner_join(tx_sizes, by = "TXNAME")
  
  # Keep longest transcript per gene
  txname_keep <- gene2txname %>%
    dplyr::group_by(GENEID) %>%
    dplyr::arrange(dplyr::desc(tx_size), .by_group = TRUE) %>%
    dplyr::summarise(TXNAME = TXNAME[[1]], .groups = "drop_last") %>%
    dplyr::pull(TXNAME)
  cds_tx <- cds_tx[names(cds_tx) %in% txname_keep]
  
  return(cds_tx)
}


#' Get the splice site locations within genes
#' 
#' This returns the position of the splice site locations,
#' within the cds of the gene.
#'
#' @param cds_tx GRangesList object containing the cds ranges of genes
#'
#' @return List containing the splice site locs in genes
#' @noRd
#'
.get_splice_genelocs = function(cds_tx){
  exon_ends <- cumsum(width(cds_tx))
  splice_start <- as.list(exon_ends - 1)
  splice_end <- as.list(exon_ends + 2)
  splice_locs_genes <- purrr::map2(splice_start, splice_end, .get_splice_single_genelocs)
  return(splice_locs_genes)
}

#' Get the splice site locations within one gene
#'
#' @param starts Start position of splice site
#' @param ends End position of splice site
#'
#' @return Vector containing the splice site locs in one gene
#' @noRd
#'
.get_splice_single_genelocs = function(starts, ends){
  splice_locs_l <- purrr::map2(starts, ends, ~seq(.x, .y))
  splice_locs <- unlist(splice_locs_l)
  return(splice_locs)
}

#' Get the potential damage per mutational context
#'
#' @param i Index of the mutational contexts
#' @param contexts_tb A tibble containing the mutational contexts
#' @param seqs DNAStringSet containing the cds sequences
#' @param verbose Boolean. Determines whether progress is printed. (Default: FALSE)
#' @param blosum62 Blosum62 matrix
#' @param splice_genelocs List containing the splice site locs in genes
#'
#' @return A tibble with the ratio of 'stop gain', 'mismatch' and 'synonymous' mutations
#' for one mutation context.
#' @noRd
#'
.single_context_damage_analysis <- function(i, contexts_tb, seqs, verbose, blosum62, splice_genelocs) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  alt_base <- context <- n <- NULL

  # Get data from this context
  contexts_tb <- contexts_tb[i, ]
  ori_bases <- contexts_tb$ori_bases
  ref_base <- contexts_tb$ref_base
  alt_bases <- contexts_tb$alt_bases[[1]]
  l_context <- contexts_tb$l_context
  r_context <- contexts_tb$r_context
  mut_pos <- contexts_tb$mut_pos
  rev_mut_pos <- contexts_tb$rev_mut_pos
  
  # Count muttypes for forward context
  muttype_counts <- .single_context_damage_analysis_strand(ori_bases, 
                                                           ref_base, 
                                                           alt_bases, 
                                                           seqs, 
                                                           mut_pos,
                                                           blosum62,
                                                           splice_genelocs) %>%
    dplyr::mutate(context = paste0(l_context, "[", ref_base, ">", alt_base, "]", r_context)) %>%
    dplyr::select(type, context, n, blosum62)

  # Get reverse context
  rev_ori_bases <- ori_bases %>%
    Biostrings::DNAString() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  rev_ref_base <- ref_base %>%
    Biostrings::DNAString() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  rev_alt_bases <- alt_bases %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>%
    as.character()

  # Count muttypes for reverse context
  muttype_counts_rev <- .single_context_damage_analysis_strand(rev_ori_bases, 
                                                               rev_ref_base, 
                                                               rev_alt_bases, 
                                                               seqs,
                                                               rev_mut_pos,
                                                               blosum62,
                                                               splice_genelocs)

  # Combine forward and reverse context
  muttype_counts$n <- muttype_counts$n + muttype_counts_rev$n
  muttype_counts$blosum62 <- muttype_counts$blosum62 + muttype_counts_rev$blosum62
  
  # Normalize
  norm_muttype_counts <- muttype_counts %>%
    dplyr::mutate(blosum62 = ifelse(type == "Missense", 
                                    blosum62 / n,
                                    NA)) %>% 
    dplyr::select(type, context, n, blosum62)

  if (verbose) {
    message(paste0("Finished with the ", ori_bases, " context."))
  }

  return(norm_muttype_counts)
}

#' Get the potential damage per mutational context for a single strand
#'
#' @param ori_bases Mutational context
#' @param ref_base Reference base
#' @param alt_bases Vector of possible alternative bases
#' @param seqs DNAStringSet containing the cds sequences
#' @param mut_pos Mutation position within context
#' @param blosum62 Blosum62 matrix
#' @param splice_genelocs List containing the splice site locs in genes
#'
#' @return A tibble with the number of 'stop gain', 'mismatch' and 'synonymous' mutations
#' for one mutation context on one strand.
#' @noRd
#'
.single_context_damage_analysis_strand <- function(ori_bases, 
                                                   ref_base, 
                                                   alt_bases, 
                                                   seqs, 
                                                   mut_pos, 
                                                   blosum62,
                                                   splice_genelocs) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  loc <- dna <- NULL

  # Determine reference codons and mutation location
  ori_bases_biostring <- Biostrings::DNAString(ori_bases)
  ref_mut_loc_l <- purrr::map2(as.list(seqs), splice_genelocs, .get_ref_codons, ori_bases_biostring, mut_pos)

  # Get ref codons from list
  ref_codons_l <- purrr::map(ref_mut_loc_l, "ref_codons")
  names(ref_codons_l) <- NULL
  ref_codons <- do.call(c, ref_codons_l)

  # Get mutation locations in codons from list
  mut_loc_in_codon_l <- purrr::map(ref_mut_loc_l, "mut_loc_in_codon")
  names(mut_loc_in_codon_l) <- NULL
  mut_loc_in_codon <- do.call(c, mut_loc_in_codon_l)


  # Count how often each combination of codon and mutated base location occurs
  tb <- tibble::tibble("loc" = mut_loc_in_codon, "dna" = as.vector(ref_codons))
  counts <- tb %>%
    dplyr::group_by(loc, dna) %>%
    dplyr::count() %>%
    dplyr::ungroup()

  # Calculate the occuring mismatch for each combi of codon and mut base location.
  muttype_counts <- purrr::map(alt_bases, .calculate_mismatches, counts, blosum62) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(ori_bases = ori_bases, ref_base = ref_base)

  return(muttype_counts)
}

#' Get reference codons for one gene
#' 
#' The reference codons are determined for one context
#' and one gene. Splice site mutations are removed.
#'
#' @param seq DNAString containing the cds for one gene
#' @param splice_single_genelocs vector containing the splice site locs in genes
#' @param ori_bases Mutational context
#' @param mut_pos Mutation position within context
#'
#' @return List. Containing reference codons and
#' the position of the possible mutation in the codon.
#' @noRd
#'
.get_ref_codons <- function(seq, splice_single_genelocs, ori_bases, mut_pos) {

  # Determine locations of context in dna
  locs <- Biostrings::matchPattern(ori_bases, seq)

  # Determine locations of mut in dna
  locs_mutbase <- start(locs) + mut_pos - 1
  
  # Remove splice site mutations
  locs_mutbase <- locs_mutbase[!locs_mutbase %in% splice_single_genelocs]
  
  # Get the reference codons
  exon_codons <- Biostrings::codons(seq)
  codon_nr <- ceiling(locs_mutbase / 3)
  ref_codons <- Biostrings::DNAStringSet(exon_codons[codon_nr])

  # Determine mutation location in reference
  mut_loc_in_codon <- dplyr::case_when(
    locs_mutbase %% 3 == 0 ~ 3,
    locs_mutbase %% 3 == 2 ~ 2,
    locs_mutbase %% 3 == 1 ~ 1
  )

  return(list("ref_codons" = ref_codons, "mut_loc_in_codon" = mut_loc_in_codon))
}

#' Calculate the possible mismatches for each of the codons
#' for one possible alternative base.
#'
#' @param alt_base Alternative base
#' @param counts Tibble of all codons.
#'
#' @return A tibble with the number of 'stop gain', 'mismatch' and 'synonymous' mutations
#' for one mutation context on one strand for one alternative base.
#' @noRd
#'
.calculate_mismatches <- function(alt_base, counts, blosum62) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  n <- ref_aa <- mut_aa <- NULL

  ref_codons_sum <- Biostrings::DNAStringSet(counts$dna)
  mut_loc_in_codon_sum <- counts$loc

  # Calculate mutated codons.
  mut_codons <- purrr::map2(as.list(ref_codons_sum), mut_loc_in_codon_sum, .mutate_codon, alt_base) %>%
    Biostrings::DNAStringSet()

  # Translate reference codons
  counts$ref_aa <- ref_codons_sum %>%
    Biostrings::translate() %>%
    as.vector()

  # Translate mutated codons
  counts$mut_aa <- mut_codons %>%
    Biostrings::translate() %>%
    as.vector()

  # Identify stop_gain, missense and synonymous mutations
  counts <- counts %>%
    dplyr::mutate(type = dplyr::case_when(
      mut_aa == "*" & ref_aa != "*" ~ "Stop_gain",
      mut_aa != ref_aa ~ "Missense",
      mut_aa == ref_aa ~ "Synonymous"
    ))
  
  # Add the BLOSUM scores
  blosum_index <- counts %>% 
    dplyr::select(ref_aa, mut_aa) %>% 
    as.matrix()
  
  counts$blosum62 <- blosum62[blosum_index] * counts$n

  # Count the number of stop_gain, missense and synonymous.
  # Also sum up the PAM scores.
  counts <- counts %>%
    dplyr::mutate(type = factor(type, levels = c("Stop_gain", "Missense", "Synonymous"))) %>%
    dplyr::group_by(type, .drop = FALSE) %>%
    dplyr::summarise(n = sum(n), blosum62 = sum(blosum62), .groups = "drop_last") %>%
    dplyr::mutate(alt_base = alt_base)
  
  
  return(counts)
}

#' Mutate codons with an alternative base
#'
#' @param ref_codon A reference codon
#' @param mut_loc_in_codon The location of the mutation in the codon
#' @param alt_base The alternative base that will be inserted
#'
#' @return A mutated version of the codon
#' @noRd
#'
.mutate_codon <- function(ref_codon, mut_loc_in_codon, alt_base) {
  ref_codon[mut_loc_in_codon] <- alt_base
  return(ref_codon)
}

#' Determine potential splice site damage
#' 
#' @param contexts_tb A tibble containing the mutational contexts
#' @param cds_tx GRangesList object containing the cds ranges of genes
#' @param ref_genome BSgenome reference genome object
#'
#' @return tibble containing splice site damage
#' @noRd
#'
.potential_splice_site_damage = function(contexts_tb, cds_tx, ref_genome){
  
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  alt_bases <- l_context <- r_context <- ref_base <- NULL
  context <- n <- blosum62 <- NULL
  
  # Get splice site sequences
  context_l <- max(c(stringr::str_length(contexts_tb$l_context), 
                     stringr::str_length(contexts_tb$r_context)))
  splice_seqs <- .get_splice_site_sequences(cds_tx, ref_genome, context_l)
  
  
  # Determine nr matched splices per context
  nr_matched_splices <- .find_potential_splice_site_muts(contexts_tb, splice_seqs)
  
  # Separate alternative bases
  splice_muts_tb <- contexts_tb %>%
    dplyr::mutate(n = nr_matched_splices) %>% 
    tidyr::unnest(alt_bases) %>% 
    dplyr::mutate(context = paste0(l_context, "[", ref_base, ">", alt_bases, "]", r_context),
                  type = "splice_site",
                  blosum62 = NA) %>% 
    dplyr::select(type, context, n, blosum62)
  
  return(splice_muts_tb)
}

#' Finds potential splice site muts based on contexts and splice site sequences
#' 
#'
#' @param contexts_tb A tibble containing the mutational contexts
#' @param splice_seqs DNAStringSet of splice site sequences
#'
#' @return Vector with the number of potential splice mutations per context
#' @noRd
#'
.find_potential_splice_site_muts = function(contexts_tb, splice_seqs){
  
  # Get unique splice site sequences
  splice_seqs_table <- BiocGenerics::table(splice_seqs)
  splice_seqs <- Biostrings::DNAStringSet(names(splice_seqs_table))
  nr_splice_seqs = as.vector(splice_seqs_table)
  
  #Forward context
  ori_bases <- Biostrings::DNAStringSet(contexts_tb$ori_bases)
  nr_matched_splices <- .find_potential_splice_site_muts_strand(ori_bases, 
                                                                splice_seqs, 
                                                                nr_splice_seqs)
  #Reverse context
  rev_ori_bases <- Biostrings::reverseComplement(ori_bases)
  nr_matched_splices_rev <- .find_potential_splice_site_muts_strand(rev_ori_bases, 
                                                                    splice_seqs, 
                                                                    nr_splice_seqs)
  #Combine both contexts
  nr_matched_splices <- nr_matched_splices + nr_matched_splices_rev
  return(nr_matched_splices)
}

#' Finds potential splice site muts based on forward/reverse context and splice site sequences
#'
#' @param ori_bases Mutational context
#' @param splice_seqs DNAStringSet of splice site sequences
#' @param nr_splice_seqs Vector with the number of potential splice mutations per context
#'
#' @return Vector with the number of potential splice mutations per context for one strand
#' @noRd
#'
.find_potential_splice_site_muts_strand = function(ori_bases, splice_seqs, nr_splice_seqs){
  nr_matched_splices <- purrr::map(as.list(ori_bases), Biostrings::vmatchPattern, splice_seqs) %>% 
    purrr::map(S4Vectors::elementNROWS) %>%
    purrr::map(function(x) x > 0) %>% 
    purrr::map_dbl(function(x) sum(x * nr_splice_seqs))
  return(nr_matched_splices)
}

#' Get splice site sequences for supplied genes
#' 
#' The 4 bases around the splice site plus the length of
#' the context is taken. This way potential mutations in 
#' the 4 bases can be found.
#'
#' @param cds_tx GRangesList object containing the cds ranges of genes
#' @param ref_genome BSgenome reference genome object
#' @param context_l Mutation context length
#'
#' @return DNAStringSet containing the splice site sequences
#' @noRd
#'
.get_splice_site_sequences = function(cds_tx, ref_genome, context_l) {
  
  # Get cds GRanges
 splice_grl = .get_splice_site_ranges(cds_tx, context_l)
  
  # Get sequences (per cds per transcript.)
  seqs <- Biostrings::getSeq(ref_genome, splice_grl) %>% 
    unlist()
  
  return(seqs)
}

#' Get splice site ranges for supplied genes.
#'
#' @param cds_tx GRangesList object containing the cds ranges of genes
#' @param context_l Mutation context length
#'
#' @return GRangesList containing the splice site ranges of genes
#' @noRd
#'
.get_splice_site_ranges = function(cds_tx, context_l){
  splice_grl <- purrr::map(as.list(cds_tx), 
                           .get_splice_site_ranges_gr,
                           context_l) %>% 
    GenomicRanges::GRangesList()
  return(splice_grl)
}

#' Get splice site ranges for one gene.
#'
#' @param cds_gr GRanges containing the cds ranges of one gene
#' @param context_l Mutation context length
#'
#' @return GRanges containing the splice site ranges of one gene
#' @noRd
#'
.get_splice_site_ranges_gr = function(cds_gr, context_l){
  
  # Create separate GRanges for start and end of splice sites.
  start_splice <- end_splice <- cds_gr
  
  # Set coordinates to the 4 bases around the splice site plus context
  BiocGenerics::end(start_splice) <- BiocGenerics::start(start_splice) + 1 + context_l
  BiocGenerics::start(start_splice) <- BiocGenerics::start(start_splice) - 2 - context_l
  BiocGenerics::start(end_splice) <- BiocGenerics::end(end_splice) - 1 - context_l
  BiocGenerics::end(end_splice) <- BiocGenerics::end(end_splice) + 2 + context_l
  
  # Remove start and end of gene, because there is no splicing there.
  start_splice <- start_splice[-1]
  end_splice <- end_splice[-length(end_splice)]
  
  #Combine start and end splice sites.
  splices_gr = c(start_splice, end_splice)
  
  return(splices_gr)
}

