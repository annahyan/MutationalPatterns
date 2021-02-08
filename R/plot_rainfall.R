#' Plot genomic rainfall
#'
#' Rainfall plot visualizes the types of mutations and intermutation distance
#' @details
#' Rainfall plots can be used to visualize the distribution of mutations
#' along the genome or a subset of chromosomes. The distance of a mutation
#' with the mutation prior to it (the intermutation distance) is plotted on
#' the y-axis on a log scale. The input GRanges are sorted before plotting.
#'
#' The colour of the points indicates the base substitution type.
#' Clusters of mutations with lower intermutation distance represent mutation
#' hotspots.
#'
#' @param vcf GRanges object
#' @param chromosomes Vector of chromosome/contig names of the reference
#' genome to be plotted
#' @param title Optional plot title
#' @param colors Vector of 6 colors used for plotting
#' @param cex Point size
#' @param cex_text Text size
#' @param ylim Maximum y value (genomic distance)
#' @return Rainfall plot
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' # Specify chromosomes of interest.
#' chromosomes <- names(genome(vcfs[[1]])[1:22])
#'
#' ## Do a rainfall plot for all chromosomes:
#' plot_rainfall(vcfs[[1]],
#'   title = names(vcfs[1]),
#'   chromosomes = chromosomes,
#'   cex = 1
#' )
#'
#' ## Or for a single chromosome (chromosome 1):
#' plot_rainfall(vcfs[[1]],
#'   title = names(vcfs[1]),
#'   chromosomes = chromosomes[1],
#'   cex = 2
#' )
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

plot_rainfall <- function(vcf,
                          chromosomes,
                          title = "",
                          colors = NA,
                          cex = 2.5,
                          cex_text = 3,
                          ylim = 1e+08) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  location <- NULL

  # If colors parameter not provided, set to default colors
  if (.is_na(colors)) {
    colors <- COLORS6
  }

  # Check color vector length
  if (length(colors) != 6) {
    stop("colors vector length not 6")
  }
  # get chromosome lengths of reference genome
  chr_length <- GenomeInfoDb::seqlengths(vcf)
  # Check for missing seqlengths
  if (sum(is.na(GenomeInfoDb::seqlengths(vcf))) > 1) {
    stop(paste(
      "Chromosome lengths missing from vcf object.\n",
      "Likely cause: contig lengths missing from the header of your vcf file(s).\n",
      "Please evaluate: seqinfo(vcf)\n",
      "To add seqlengths to your vcf GRanges object use: seqlengths(vcf) <- my_chr_lengths" 
    ), call. = FALSE)
  }

  # Sort the input
  vcf <- BiocGenerics::sort(vcf)

  # subset on the chromosomes selected by the user.
  chr_length <- chr_length[names(chr_length) %in% chromosomes]

  # cumulative sum of chromosome lengths
  chr_cum <- c(0, cumsum(as.numeric(chr_length)))

  # Plot chromosome labels without "chr"
  names(chr_cum) <- names(chr_length)
  labels <- gsub("chr|chromosome_|chromosome|group|group_|chrom", 
                 "", 
                 names(chr_length), 
                 ignore.case = TRUE)

  # position of chromosome labels.
  # Calculated by taking the average between two adjacent chr_cums.
  m <- (chr_cum[-1] + chr_cum[-length(chr_cum)]) / 2

  # Determine mutation characteristics per chromosome.
  tb_l <- purrr::map(chromosomes, function(chr) {

    # Subset variants to chromosome
    chr_subset <- vcf[GenomeInfoDb::seqnames(vcf) == chr]

    # If there are not enough variants, then an empty tibble is returned.
    n <- length(chr_subset)
    if (n <= 1) {
      tb <- tibble::tibble(
        "type" = character(0),
        "location" = double(0),
        "distance" = integer(0),
        "chromosome" = character(0)
      )
      return(tb)
    }

    # Determine type, location and distance to previous mut.
    type <- mut_type(chr_subset)[-1]
    loc <- (start(chr_subset) + chr_cum[chr])[-1]
    dist <- diff(start(chr_subset))

    # Combine all into a tibble.
    tb <- tibble::tibble(
      "type" = type,
      "location" = loc,
      "distance" = dist,
      "chromosome" = chr
    )
    return(tb)
  })

  # Combine the different chromosomes.
  data <- do.call(rbind, tb_l)

  # Removes colors based on missing mutation types.  This prevents colors from
  # shifting when comparing samples with low mutation counts.
  typesin <- SUBSTITUTIONS %in% unique(data$type)
  colors <- colors[typesin]

  # make rainfall plot
  plot <- ggplot(data, aes(x = location, y = distance)) +
    geom_point(aes(colour = factor(type)), cex = cex) +
    geom_vline(xintercept = as.vector(chr_cum), linetype = "dotted") +
    annotate("text", x = m, y = ylim, label = labels, cex = cex_text) +
    scale_y_log10() +
    scale_colour_manual(values = colors) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(chr_cum))) +
    labs(x = "Genomic Location", y = "Genomic Distance", title = title) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    guides(colour = guide_legend(nrow = 1))

  return(plot)
}
