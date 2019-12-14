#' Plot signature contribution barplot
#' 
#' Plot contribution of signatures
#' 
#' @param contribution Signature contribution matrix
#' @param signatures Signature matrix
#' @param index optional sample subset parameter
#' @param coord_flip Flip X and Y coordinates, default = FALSE
#' @param mode "relative" or "absolute"; to plot the relative contribution or
#' absolute number of mutations, default = "relative"
#' @param palette A color palette like c("#FF0000", "#00FF00", "9999CC") that
#' will be used as colors in the plot.  By default, ggplot2's colors are used
#' to generate a palette.
#'
#' @return Stacked barplot with contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Optionally set column and row names.
#' colnames(nmf_res$signatures) = c("Signature A", "Signature B")
#' rownames(nmf_res$contribution) = c("Signature A", "Signature B")
#'
#' ## The following are examples of contribution plots.
#' plot_contribution(nmf_res$contribution,
#'                     nmf_res$signature,
#'                     mode = "relative")
#' 
#' plot_contribution(nmf_res$contribution,
#'                     nmf_res$signature,
#'                     mode = "absolute")
#' 
#' plot_contribution(nmf_res$contribution,
#'                     nmf_res$signature,
#'                     mode = "absolute",
#'                     index = c(1,2))
#' 
#' plot_contribution(nmf_res$contribution,
#'                     nmf_res$signature,
#'                     mode = "absolute",
#'                     coord_flip = TRUE)
#'
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}}
#'
#' @export

plot_contribution = function(contribution,
                                signatures,
                                index=c(),
                                coord_flip=FALSE,
                                mode="relative",
                                palette=c())
{
    # check mode parameter
    if(!(mode == "relative" | mode == "absolute"))
        stop("mode parameter should be either 'relative' or 'absolute'")

    # optional subsetting if index parameter is provided
    if(length(index > 0)){contribution = contribution[,index]}

    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    Sample = NULL
    Contribution = NULL
    Signature = NULL

    if (mode == "relative")
    {
        # Plot contribution
        m_contribution = melt(contribution)
        colnames(m_contribution) = c("Signature", "Sample", "Contribution")

        plot = ggplot(m_contribution,
                        aes(x = factor(Sample),
                            y = Contribution,
                            fill = factor(Signature),
                            order = Sample)) +
            geom_bar(position = "fill", stat="identity", colour="black")  +
            # ylabel
            labs(x = "", y = "Relative contribution") +
            # white background
            theme_bw() +
            # no gridlines
            theme(panel.grid.minor.x=element_blank(),
                    panel.grid.major.x=element_blank()) +
            theme(panel.grid.minor.y=element_blank(),
                    panel.grid.major.y=element_blank())
    }

    # Handle the absolute mode.
    else 
    {
        if(missing(signatures))
            stop(paste("For contribution plotting in mode 'absolute':",
                        "also provide signatures matrix"))

        # total number of mutations per siganture
        total_signatures = colSums(signatures) 

        # calculate signature contribution in absolute number of signatures
        abs_contribution = contribution * total_signatures

        # Plot contribution
        m_contribution = melt(abs_contribution)
        colnames(m_contribution) = c("Signature", "Sample", "Contribution")

        plot = ggplot(m_contribution, aes(x = factor(Sample),
                                            y = Contribution,
                                            fill = factor(Signature),
                                            order = Sample)) + 
            geom_bar(stat="identity", colour = "black")  +  
            # ylabel
            labs(x = "", y = "Absolute contribution \n (no. mutations)") +  
            # white background
            theme_bw() +
            # no gridlines
            theme(panel.grid.minor.x=element_blank(),
                    panel.grid.major.x=element_blank()) +
            theme(panel.grid.minor.y=element_blank(),
                    panel.grid.major.y=element_blank())
    }

    # Allow custom color palettes.
    if (length(palette) > 0)
        plot = plot + scale_fill_manual(name="Signature", values=palette)
    else
        plot = plot + scale_fill_discrete(name="Signature")

    # Handle coord_flip.
    if (coord_flip)
        plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
    else
        plot = plot + xlim(levels(factor(m_contribution$Sample)))
                
    return(plot)
}
