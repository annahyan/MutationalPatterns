#' Plot a riverplot
#' 
#' Function to plot a SNV mutation matrix as a riverplot.
#' This is especially useful when looking at a wide
#' mutational context
#'
#' @param mut_matrix Matrix containing mutation counts.
#' @param condensed More condensed plotting format. Default = F.
#' 
#' @return A ggplot object
#' 
#' @export
#' @import ggplot2
#' @import ggalluvial
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{plot_96_profile}},
#' \code{\link{plot_profile_heatmap}}
#' 
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
#' plot_river(mut_mat[,c(1,4)])
#'
#' ## Get extended matrix
#' mut_mat_extended <- readRDS(system.file("states/mut_mat_data_extended.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create heatmap of extended profile
#' plot_river(mut_mat_extended[,c(1,4)])
#' 
#' ## Create condensed version of riverplot
#' plot_river(mut_mat_extended[,c(1,4)], condensed = TRUE)
#' 
plot_river = function(mut_matrix, condensed = FALSE){
    
    context_m <- .river_create_context_matrix(mut_matrix)
    lodes_tb <- .river_create_lodes(mut_matrix, context_m)
    fig <- .river_plot_river(lodes_tb, mut_matrix, condensed)
    
    return(fig)
}

#' Split the contexts of a mutation matrix
#'
#' Function to split the contexts of a mutation matrix
#' in multiple columns.
#' 
#' @param mut_matrix Matrix containing mutation counts.
#'
#' @return Matrix containing the contexts of a mutation matrix.
#' @noRd
#' @importFrom magrittr %>%
#' 
.river_create_context_matrix = function(mut_matrix){
    
    # Split left and right context from mutations
    context_m <- mut_matrix %>% 
        rownames() %>% 
        stringr::str_split("\\[|\\]", simplify = TRUE)
    
    # Split contexts into single bases
    left_m <- stringr::str_split(context_m[,1], "", simplify = TRUE)
    right_m <- stringr::str_split(context_m[,3], "", simplify = TRUE)
    
    # Combine all
    context_m <- cbind(left_m, context_m[,2], right_m)
    
    # Set column names
    mut_positions <- c((-1*rev(seq_len(ncol(left_m)))), 0, seq_len(ncol(right_m)))
    colnames(context_m) <- paste0("pos_", mut_positions)
    
    return(context_m)
}


#' Create a long format lodes tibble.
#'
#' @param mut_matrix Matrix containing mutation counts.
#' @param context_m Matrix containing the contexts of a mutation matrix.
#'
#' @return A tibble in lodes format containing the mutation contexts.
#' @noRd
#'
#' @importFrom magrittr %>%
#' 
.river_create_lodes = function(mut_matrix, context_m){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    pos <- type <- sample_name <- rowid <- nrmuts <- NULL
    
    mut_df <- as.data.frame(mut_matrix)
    context_df <- as.data.frame(context_m)
    
    # Transform data into long (lodes) format.
    lodes_tb <- cbind(mut_df, context_df) %>% 
        tibble::rowid_to_column() %>% 
        tidyr::pivot_longer(cols = dplyr::contains("pos_"), names_to = "pos", values_to = "type") %>% 
        tidyr::pivot_longer(cols = c(-pos, -type, -rowid), 
                            names_to = "sample_name", 
                            values_to = "nrmuts") %>% 
        dplyr::mutate(pos = stringr::str_remove(pos, "pos_"),
            pos = factor(pos, levels = unique(pos)),
            type = factor(type, levels = unique(type)),
            sample_name = factor(sample_name, levels = unique(sample_name)))
    return(lodes_tb)
}

#' Plot a riverplot with a lodes tibble as input
#'
#' @param lodes_tb A tibble in lodes format containing the mutation contexts.
#' @param condensed More condensed plotting format.
#'
#' @return A ggplot object
#' @noRd
#'
#' @import ggplot2
#' @import ggalluvial
#' 
.river_plot_river = function(lodes_tb, mut_matrix, condensed){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    rowid <- nrmuts <- NULL
    
    
    # Set colours
     colours_v <- c("#EF2D56", "#446DF6","#6FAB94", "#FFBC63", 
                  "#2EBAED", "#000000", "#DE1C14", "#D4D2D2", 
                  "#ADCC54", "#F0D0CE")
    types <- c("T", "C", "A", "G", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    typesin <- types %in% unique(lodes_tb$type)
    used_colours <- colours_v[typesin]
    
    # Change plotting parameters based on whether plot should be condensed.
    if (condensed == TRUE) {
        spacing <- 0
    } else {
        spacing <- 0.5
    }
    
    # Create facet texts
    nr_muts <- colSums(mut_matrix)
    facet_labs_y <- stringr::str_c(colnames(mut_matrix), " (n = ", nr_muts, ")")
    names(facet_labs_y) <- colnames(mut_matrix)
    
    #Add stratum stat.
    StatStratum <- ggalluvial::StatStratum
    
    # Create plot 
    fig <- ggplot(lodes_tb, aes(x = pos, 
                  stratum = type, 
                  y = nrmuts, 
                  alluvium = rowid, 
                  fill = type, 
                  label = type)) + 
        geom_stratum() + 
        geom_flow() + 
        geom_text(stat = StatStratum, size = 3, colour = "white") +
        facet_grid(sample_name ~ ., scales = "free_y",
                   labeller = labeller(sample_name = facet_labs_y)) +
        scale_fill_manual(values = used_colours) +
        labs(x = "Position", y = "Nr mutations") +
        theme_bw() +
        theme(panel.spacing.y = unit(spacing, "lines"))
    
    return(fig)
}
