#' Plot strand per base substitution type
#' 
#' For each base substitution type and transcriptional strand the total number
#' of mutations and the relative contribution within a group is returned.
#' @param strand_bias_df data.frame, result from strand_bias function
#' @param mode Either "absolute" for absolute number of mutations, or
#' "relative" for relative contribution, default = "relative"
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' strand_counts = strand_occurrences(mut_mat_s, by=tissue)
#'
#' #' ## Plot the strand in relative mode.
#' strand_plot = plot_strand(strand_counts)
#'
#' #' ## Or absolute mode.
#' strand_plot = plot_strand(strand_counts, mode = "absolute")
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

plot_strand = function(strand_bias_df, mode = "relative", colors)
{
    # if colors parameter not provided, set to default colors
    if (missing(colors))
        colors=COLORS6

    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    type = NULL
    relative_contribution = NULL
    no_mutations = NULL

    # Plot relative contribution within each group
    if(mode == "relative")
    {
        plot = ggplot(strand_bias_df, aes(x=type,
                                            y=relative_contribution,
                                            fill=type,
                                            alpha=strand)) +
            geom_bar(stat="identity",
                        position = "dodge",
                        colour="black",
                        cex=0.5) + 
            scale_fill_manual(values= colors) +
            scale_alpha_discrete(range = c(1, 0.4)) +
            ylab("Relative contribution") +
            facet_grid(. ~ group) +
            theme_bw() +
            scale_x_discrete(breaks=NULL) +
            xlab("")
    }

    # Plot absolute contribution within each group
    else if (mode == "absolute")
    {
        plot = ggplot(strand_bias_df, aes(x=type,
                                            y=no_mutations,
                                            fill=type,
                                            alpha=strand)) +
            geom_bar(stat="identity",
                        position = "dodge",
                        colour="black",
                        cex=0.5) +
            scale_fill_manual(values= colors) +
            scale_alpha_discrete(range = c(1, 0.4)) +
            ylab("Total number of mutations") +
            facet_grid(. ~ group) +
            theme_bw() +
            scale_x_discrete(breaks=NULL) +
            xlab("")
    }

    return(plot)
}
