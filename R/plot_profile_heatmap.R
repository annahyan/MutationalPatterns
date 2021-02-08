#' Plot a mutation matrix as a heatmap
#'
#' Function to plot a SNV mutation matrix as a heatmap.
#' This is especially useful when looking at a wide mutational context.
#'
#'
#' @param mut_matrix Matrix containing mutation counts.
#' @param by Optional grouping variable
#' @param max Maximum value used for plotting the relative contributions.
#'   Contributions that are higher will have the maximum colour. (Default: 0.02)
#' @param condensed More condensed plotting format. Default = F.
#' 
#' @return A ggplot object
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{plot_96_profile}},
#' \code{\link{plot_river}}
#' @examples
#'
#' ## See the 'mut_matrix()' examples for how we obtained the
#' ## mutation matrix information:
#' ## Get regular matrix
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create heatmap of profile
#' plot_profile_heatmap(mut_mat, max = 0.1)
#'
#' ## Get extended matrix
#' mut_mat_extended <- readRDS(system.file("states/mut_mat_data_extended.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create heatmap of extended profile
#' plot_profile_heatmap(mut_mat_extended)
#'
#' ## Or plot heatmap per tissue
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' plot_profile_heatmap(mut_mat_extended, by = tissue)
#'
#' ## Or plot the heatmap per sample.
#' plot_profile_heatmap(mut_mat_extended,
#'   by = colnames(mut_mat_extended),
#'   max = 0.05
#' )
#' 
#' 
#' ## Create a condensed heatmap of extended profile
#' plot_profile_heatmap(mut_mat_extended, condensed = TRUE)

plot_profile_heatmap <- function(mut_matrix, 
                                 by = NA, 
                                 max = 0.02,
                                 condensed = FALSE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  fullcontext <- l_context <- r_context <- muttype <- NULL
  nrmuts <- rel_nrmuts <- NULL

  # check arguments
  if (!inherits(mut_matrix, "matrix")) {
    stop("mut_matrix must be a matrix", call. = FALSE)
  }
  # matrix should have row and colnames
  if (length(colnames(mut_matrix)) == 0) {
    stop("mut_matrix is missing colnames", call. = FALSE)
  }
  if (length(rownames(mut_matrix)) == 0) {
    stop("mut_matrix is missing rownames", call. = FALSE)
  }


  # Transform the data into long format and get left/right context seperately.
  tb <- mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("fullcontext") %>%
    tidyr::pivot_longer(-fullcontext, names_to = "sample", values_to = "nrmuts") %>%
    tidyr::separate("fullcontext", # Separate context into left context, mut and right context
      into = c("l_context", "muttype", "r_context"),
      sep = "\\[|\\]"
    ) %>%
    dplyr::mutate(
      mut = factor(muttype, levels = unique(muttype)),
      r_context = factor(r_context, levels = unique(r_context)),
      l_context = factor(l_context, levels = Biostrings::reverse(unique(l_context)))
    )

  # Make data relative
  tb <- tb %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(rel_nrmuts = nrmuts / sum(nrmuts)) %>% 
    dplyr::ungroup()

  # Add sample grouping
  if (.is_na(by)) {
    by <- "all"
  }
  tb_by <- tibble::tibble(
    "sample" = unique(tb$sample),
    "by" = by
  )
  tb <- tb %>%
    dplyr::left_join(tb_by, by = "sample")

  # Combine samples based on grouping
  tb <- tb %>%
    dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, l_context, muttype, r_context) %>%
    dplyr::summarise(
      nrmuts = sum(nrmuts),
      rel_nrmuts = mean(rel_nrmuts),
      total_indv = dplyr::n(),
      .groups = "drop_last"
    ) %>%
    dplyr::ungroup()

  # If value is higher than y_max, change it to ymax. (Prevents plotting issues)
  tb <- tb %>%
    dplyr::mutate(rel_nrmuts = ifelse(rel_nrmuts > max, max, rel_nrmuts))


  # Count number muts per sample_group
  tot_muts_samplegroups_tb <- tb %>%
    dplyr::group_by(by) %>%
    dplyr::summarise(nrmuts = sum(nrmuts), .groups = "drop_last")

  # Create y facet labels
  facet_labs_y <- stringr::str_c(tot_muts_samplegroups_tb$by, " (n = ", tot_muts_samplegroups_tb$nrmuts, ")")
  names(facet_labs_y) <- tot_muts_samplegroups_tb$by

  # Count number muts per mut
  tot_muts_muttype_tb <- tb %>%
    dplyr::group_by(muttype) %>%
    dplyr::summarise(nrmuts = sum(nrmuts), .groups = "drop_last")

  # Create x facet labels
  facet_labs_x <- stringr::str_c(tot_muts_muttype_tb$muttype, " (n = ", tot_muts_muttype_tb$nrmuts, ")")
  names(facet_labs_x) <- tot_muts_muttype_tb$muttype

  # Set plotting parameters
  context_size <- stringr::str_length(tb$l_context[1])
  if (context_size == 1) {
    axis_size <- 10
  } else if (context_size == 2) {
    axis_size <- 8
  } else if (context_size == 3) {
    axis_size <- 4
  } else {
    axis_size <- 3
  }

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    spacing <- 0
  } else {
    spacing <- 0.5
  }
  
  
  
  # Create plot
  fig <- ggplot(tb, aes(x = r_context, 
                        y = l_context, 
                        fill = rel_nrmuts)) +
    geom_raster() +
    scale_fill_distiller(
      palette = "YlGnBu",
      direction = 1,
      name = "Relative contribution",
      limits = c(0, max)
    ) +
    facet_grid(by ~ muttype, labeller = labeller(by = facet_labs_y, muttype = facet_labs_x)) +
    labs(x = "3' context", y = "5' context") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = axis_size),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x = unit(spacing, "lines"),
      panel.spacing.y = unit(spacing, "lines")
    )

  return(fig)
}
