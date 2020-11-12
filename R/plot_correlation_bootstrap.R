#' Plots the correlation between bootstrapped signature contributions
#'
#' This function plots the pearson correlation between signatures.
#' This can be done per sample or for all samples together.
#' It returns a list of the created figures.
#'
#' @param contri_boots A dataframe with bootstrapped signature contributions.
#' @param per_sample Whether or not a plot should be made per sample. Default: TRUE.
#'
#' @return A list of ggplot2 objects if run per sample.
#' Else it returns a single ggplot2 object.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ## Get a dataframe with bootstrapped signature contributions.
#' ## See 'fit_to_signatures_bootstrapped()' for how to do this.
#' contri_boots <- readRDS(system.file("states/bootstrapped_snv_refit.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Plot the correlations between signatures per sample
#' fig_l <- plot_correlation_bootstrap(contri_boots)
#'
#' ## Look at the figure of the first sample.
#' fig_l[[1]]
#'
#' ## You can also look at the correlation for all samples combined
#' plot_correlation_bootstrap(contri_boots, per_sample = FALSE)
plot_correlation_bootstrap <- function(contri_boots, per_sample = TRUE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  if (per_sample) {
    # Determine samples
    samples <- contri_boots %>%
      rownames() %>%
      gsub("_[^_]+$", "", .)
    unique_samples <- unique(samples)
    n_samples <- length(unique_samples)

    # Create a figure per sample
    figs <- vector("list", n_samples)
    for (i in seq_len(n_samples)) {

      # Get contributions from sample
      sample <- unique_samples[i]
      index <- sample == samples
      contri_boots_sample <- contri_boots[index, , drop = FALSE]

      # Create figure
      fig <- .plot_correlation_bootstrap_sample(contri_boots_sample, sample)
      figs[[i]] <- fig
    }
  } else {
    # Create figure for all samples combined.
    figs <- .plot_correlation_bootstrap_sample(contri_boots, "all")
  }
  return(figs)
}


#' Plots the correlation between bootstrapped signature contributions for one sample.
#'
#' @param contri_boots A dataframe with bootstrapped signature contributions.
#' @param sample The name of the sample
#'
#' @return A ggplot2 object
#'
#' @importFrom magrittr %>%
#' @noRd
#'
.plot_correlation_bootstrap_sample <- function(contri_boots, sample) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  sig_row <- sig_column <- NULL

  # Get correlations
  withCallingHandlers(
    {
      sig_cor <- cor(contri_boots)
    },
    warning = function(w) {
      if (grepl("the standard deviation is zero", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )


  # Make data long
  sig_cor_tb <- sig_cor %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sig_row") %>%
    tidyr::pivot_longer(-sig_row, names_to = "sig_column", values_to = "cor") %>%
    dplyr::mutate(
      sig_row = factor(sig_row, levels = unique(sig_row)),
      sig_column = factor(sig_column, levels = unique(sig_column))
    )

  # Create figure
  fig <- ggplot(data = sig_cor_tb, aes(x = sig_column, y = sig_row, fill = cor), order = NULL) +
    geom_raster() +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Correlation",
                         na.value = "grey85") +
    labs(x = NULL, y = NULL, title = sample) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      text = element_text(size = 12),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  return(fig)
}
