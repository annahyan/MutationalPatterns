#' Plot 192 trinucleotide profile
#'
#' Plot relative contribution of 192 trinucleotides      
#' @param mut_matrix 192 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors 6 value color vector
#' @param condensed More condensed plotting format. Default = F.
#' @return 192 trinucleotide profile plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## mutation matrix with transcriptional strand information:
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Extract the signatures.
#' ## This is a computationally intensive task, so we load a precomputed
#' ## version instead.
#' # nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
#' nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Optionally, provide signature names
#' colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")
#'
#' ## Generate the plot
#' plot_192_profile(nmf_res_strand$signatures)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{extract_signatures}}
#'
#' @export

plot_192_profile = function(mut_matrix, colors, ymax = 0.2, condensed = FALSE)
{
    # Relative contribution
    norm_mut_matrix = apply(mut_matrix, 2, function(x) x / sum(x))

    # Check color vector length
    # Colors for plotting
    if(missing(colors)){colors=COLORS6}
    if(length(colors) != 6){stop("Provide colors vector with length 6")}
    context = rep(CONTEXTS_96, each=2)
    substitution = rep(SUBSTITUTIONS, each=32)
    # get strand from rownames of mut_matrix
    strand = sapply(rownames(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
    
    # Replace mutated base with dot to get context
    substring(context, 2, 2) = "."
    
    # Construct dataframe
    df = data.frame(substitution = substitution,
                    context = context,
                    strand = strand)

    rownames(norm_mut_matrix) = NULL

    df2 = cbind(df, as.data.frame(norm_mut_matrix))
    df3 = melt(df2, id.vars = c("substitution", "context", "strand"))

    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    value = NULL
    
    if (condensed)
    {
      plot = ggplot(data=df3, aes(x=context,
                                  y=value,
                                  fill=substitution,
                                  width=1,
                                  alpha=strand)) +
        geom_bar(stat="identity", colour="black", size=.2) +
        scale_fill_manual(values=colors) +
        facet_grid(variable ~ substitution) +
        ylab("Relative contribution") +
        coord_cartesian(ylim=c(0,ymax)) +
        scale_y_continuous(breaks=seq(0, ymax, 0.1)) +
        # no legend
        guides(fill=FALSE) +
        # white background
        theme_bw() +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
              axis.text.y=element_text(size=8),
              axis.title.x=element_text(size=12),
              axis.text.x=element_text(size=5,angle=90,vjust=0.4),
              strip.text.x=element_text(size=9),
              strip.text.y=element_text(size=9),
              panel.grid.major.x = element_blank(),
              panel.spacing.x = unit(0, "lines"))
    } else {
        plot = ggplot(data=df3, aes(x=context,
                                y=value,
                                fill=substitution,
                                width=0.6,
                                alpha=strand)) +
        geom_bar(stat="identity", colour="black", size=.2) + 
        scale_fill_manual(values=colors) + 
        facet_grid(variable ~ substitution) + 
        ylab("Relative contribution") + 
        coord_cartesian(ylim=c(0,ymax)) +
        scale_y_continuous(breaks=seq(0, ymax, 0.1)) +
        # no legend
        guides(fill=FALSE) + 
        # white background
        theme_bw() +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=12),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=9),
                strip.text.y=element_text(size=9),
                panel.grid.major.x = element_blank())
    }
    
    return(plot)
}
