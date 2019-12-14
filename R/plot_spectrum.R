#' Plot point mutation spectrum
#'    
#' @param type_occurrences Type occurrences matrix
#' @param CT Distinction between C>T at CpG and C>T at other
#' sites, default = FALSE
#' @param by Optional grouping variable
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = TRUE
#' @return Spectrum plot
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#' @importFrom plyr ddply
#' @importFrom plyr summarise
#' @importFrom stats sd
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#' 
#'
#' ## Load a reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences = mut_type_occurrences(vcfs, ref_genome)
#' 
#' ## Plot the point mutation spectrum over all samples
#' plot_spectrum(type_occurrences)
#'
#' ## Or with distinction of C>T at CpG sites
#' plot_spectrum(type_occurrences, CT = TRUE)
#'
#' ## Or without legend
#' plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
#'
#' ## Or plot spectrum per tissue
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#;
#' plot_spectrum(type_occurrences, by = tissue, CT = TRUE)
#'
#' ## You can also set custom colors.
#' my_colors = c("pink", "orange", "blue", "lightblue",
#'                 "green", "red", "purple")
#'
#' ## And use them in a plot.
#' plot_spectrum(type_occurrences,
#'                 CT = TRUE,
#'                 legend = TRUE,
#'                 colors = my_colors)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_type_occurrences}}
#'
#' @export

plot_spectrum = function(type_occurrences, CT=FALSE, by, colors, legend=TRUE)
{
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    value = NULL
    nmuts = NULL
    sub_type = NULL
    variable = NULL
    error_pos = NULL
    stdev = NULL

    # If colors parameter not provided, set to default colors
    if (missing(colors))
        colors = COLORS7

    # Check color vector length
    if (length(colors) != 7)
        stop("Colors parameter: supply color vector with length 7")

    # Distinction between C>T at CpG or not
    if (CT == FALSE)
        type_occurrences = type_occurrences[,1:6] 
    else
        type_occurrences = type_occurrences[,c(1:2,8,7,4:6)]

    # Relative contribution per sample
    df2 = type_occurrences / rowSums(type_occurrences)

    # If grouping variable not provided, set to "all"
    if (missing(by))
        by="all"

    # Add by info to df
    df2$by = by

    # Reshape
    df3 = melt(df2, id.vars = "by")

    # Count number of mutations per mutation type
    counts = melt(type_occurrences, measure.vars = colnames(type_occurrences))
    df4 = cbind(df3, counts$value)
    colnames(df4)[4] = "nmuts" 

    # Calculate the mean and the stdev on the value for each group broken
    # down by type + variable
    x = ddply(df4, c("by", "variable"), summarise,
                mean = mean(value), stdev = sd(value))

    info_x = ddply(df4, c("by"), summarise, total_individuals = sum(value),
                    total_mutations = sum(nmuts))

    x = merge(x, info_x)
    info_type = data.frame(sub_type = c("C>A", "C>G", "C>T", "C>T",
                                        "C>T", "T>A", "T>C", "T>G"),
                            variable = c("C>A", "C>G", "C>T", "C>T at CpG",
                                        "C>T other", "T>A", "T>C", "T>G"))
    x = merge(x,info_type)
    x$total_mutations = prettyNum(x$total_mutations, big.mark = ",")
    x$total_mutations = paste("No. mutations =",
                                as.character(x$total_mutations))

    # Define positioning of error bars
    x$error_pos = x$mean

    # Define colors for plotting
    if(CT == FALSE)
        colors = colors[c(1,2,3,5:7)]

    # If C>T stacked bar (distinction between CpG sites and other)
    else
    {
        # Adjust positioning of error bars for stacked bars
        # mean of C>T at CpG should be plus the mean of C>T at other
        x = x[order(x$by),]
        CpG = which(x$variable == "C>T at CpG")
        other = which(x$variable == "C>T other")
        x$error_pos[CpG] = x$error_pos[other] + x$error_pos[CpG]

        # Define order of bars
        order = order(factor(x$variable, levels = c("C>A", "C>G", "C>T other",
                                                    "C>T at CpG","T>A", "T>C",
                                                    "T>G")))
        x = x[order,]
    }

    # Make barplot
    plot = ggplot(data=x, aes(x=sub_type,
                                y=mean,
                                fill=variable,
                                group=sub_type)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=colors, name="Point mutation type") +
        theme_bw() +
        xlab("") +
        ylab("Relative contribution") +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                panel.grid.major.x = element_blank())
    
    # check if standard deviation error bars can be plotted
    if(sum(is.na(x$stdev)) > 0)
      warning("No standard deviation error bars can be plotted, because there is only one sample per mutation spectrum")
    else
      plot = plot + geom_errorbar(aes(ymin=error_pos-stdev, 
                                      ymax=error_pos+stdev), width=0.2)
    
      
    # Facetting
    if (length(by) == 1)
        plot = plot + facet_wrap( ~ total_mutations)
    else
        plot = plot + facet_wrap(by ~ total_mutations)

    # Legend
    if (legend == FALSE)
        plot = plot + theme(legend.position="none")

    return(plot)
} 
