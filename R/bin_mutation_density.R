#' Bin the genome based on mutation density
#'
#' This function splits the genome based on the mutation density.
#' The density is calculated per chromosome. The density is split
#' into bins. The difference in density between subsequent bins is the same
#' for all bins. In other words, the difference in density between bins 1 and
#' 2 is the same as between bins 2 and 3.
#' The function returns a GRangesList. Each GRanges in the list contains the
#' regions associated with that bin. This can be used with the
#' 'split_muts_region()' function.
#'
#'
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSgenome reference genome object
#' @param nrbins The number of bins in which to separate the genome
#' @param man_dens_cutoffs Manual density cutoffs to use.
#'
#' @family genomic_regions
#' @return GRangesList
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ### See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Determine region density
#' dens_grl <- bin_mutation_density(grl, ref_genome, nrbins = 3)
#' names(dens_grl) <- c("Low", "Medium", "High")
#'
#'
#' ## You can also use manual cutoffs. This feature is meant for more
#' ## advanced users. It can be usefull if you want to find highly mutated regions, with
#' ## a consistent cutoff between analyses.
#' dens_grl_man <- bin_mutation_density(grl, ref_genome, man_dens_cutoffs = c(0, 2e-08, 1))
bin_mutation_density <- function(vcf_list, 
                                 ref_genome, 
                                 nrbins = 3,
                                 man_dens_cutoffs = NA) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }

  # Merge if grl.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    gr <- unlist(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    gr <- vcf_list
  } else {
    .not_gr_or_grl(vcf_list)
  }

  # Get ref genome
  ref_genome <- BSgenome::getBSgenome(ref_genome)

  # Sort the input
  gr <- BiocGenerics::sort(gr)

  # Determine density per chromosome
  chroms <- GenomeInfoDb::seqlevelsInUse(gr)
  dens_gr <- purrr::map(chroms, 
                        .get_mutation_density_chrom, 
                        ref_genome, 
                        gr, 
                        chroms) %>%
    do.call(c, .)

  # Set density break type
  if (!.is_na(man_dens_cutoffs)) {
    breaks <- man_dens_cutoffs
  } else {
    breaks <- nrbins
  }

  # Split density in classes
  dens_gr$dens <- cut(dens_gr$dens, breaks = breaks, labels = FALSE)
  dens_grl <- split(dens_gr, dens_gr$dens)
  dens_grl <- GenomicRanges::reduce(dens_grl)

  return(dens_grl)
}

#' Get mutation density for a single chromosome
#'
#' This function determines the mutation density for a single chromosome.
#' It returns a GRanges object, containing the densities.
#'
#' @param chrom String describing a chromosome
#' @param ref_genome BSgenome reference genome object
#' @param gr GRanges object
#' @param chroms String vector of chromosomes
#'
#' @return GRanges
#' @noRd
#'
.get_mutation_density_chrom <- function(chrom, ref_genome, gr, chroms) {

  # Select muts in chrom
  gr <- gr[GenomeInfoDb::seqnames(gr) == chrom]

  # Determine position
  half_width <- (BiocGenerics::end(gr) - BiocGenerics::start(gr)) / 2
  pos <- BiocGenerics::start(gr) + half_width

  # Calculate density. Only calculate within chromosome size
  chr_size <- GenomeInfoDb::seqlengths(ref_genome)[chrom]
  dens <- stats::density(pos, bw = "SJ", from = 1, to = chr_size)

  # Determine location of density bins.
  # The positions of the density estimates x are taken as the middle of the bins
  half_dist <- (dens$x[2] - dens$x[1]) / 2
  start <- ceiling(dens$x - half_dist)
  end <- floor(dens$x + half_dist)

  # Clip again to chromosome size
  start[1] <- 1
  end[length(end)] <- chr_size

  # Transform to GRanges
  dens_gr <- GenomicRanges::GRanges(seqnames = chrom, 
                                    ranges = IRanges::IRanges(start, end))
  dens_gr$dens <- dens$y
  GenomeInfoDb::seqlevels(dens_gr) <- chroms

  return(dens_gr)
}
