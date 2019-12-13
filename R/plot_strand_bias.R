#' Plot strand bias per base substitution type per group
#'
#' @param strand_bias A data.frame, result from strand_bias function
#' @param colors (Optional) Named list with 6 value color vector "snv" for snv, 
#' 10 value color vector "dbs" for dbs.
#' For indels give same number of colors as there are classes.
#' @param max_yaxis (Optional) Binary option to set y axis to maximum value of all
#' mutation types. Default is 'max = FALSE'
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
#' ## Perform the strand bias test.
#' strand_counts = strand_occurrences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' ## Plot the strand bias.
#' plot_strand_bias(strand_bias)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{strand_bias_test}}
#' \code{\link{plot_strand}}
#'
#' @export

plot_strand_bias = function(strand_bias, colors, max_yaxis = FALSE)
{
  if (class(strand_bias) == "list")
  { 
    strand_bias = do.call(rbind, strand_bias) 
    method = "split"
  } else { method = "combine" }
  
  # Set default colors
  if (missing(colors))
  {
    colors = c()
    if (any(grepl("snv", strand_bias$mutation))) { colors = c(colors, COLORS6) }
    if (any(grepl("dbs", strand_bias$mutation))) { colors = c(colors, COLORS10) }
      if (any(grepl("indel", strand_bias$mutation))) 
        colors = c(colors, COLORS_INDEL)
  }
  
  # get variable names
  var_names = colnames(strand_bias)[4:5]
  
  # determine max y value for plotting
  # = log2 ratio with pseudo counts of 0.1
  log2_ratio = log2(  (strand_bias[,4]+0.1) /
                        (strand_bias[,5]+0.1))
  
  # max yvalue for plotting plus
  max = round(max(abs(log2_ratio)), digits = 1) + 0.1
  pos_stars = abs(log2((strand_bias[, 4])/(strand_bias[, 5] + 0.1)))
  max_pos_star = round(max(pos_stars[is.finite(pos_stars)]), digits = 1) + 0.1
  if(max < max_pos_star){
    max = max_pos_star
  }
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL
  type = NULL
  significant = NULL
  
  # add label for infinite values
  label2 = log2(strand_bias$ratio)
  select = which(is.finite(label2))
  label2[select] = " "
  
  strand_bias$mutation = factor(strand_bias$mutation, levels = unique(strand_bias$mutation))
  strand_bias$type = factor(strand_bias$type, levels = unique(strand_bias$type))
  
  # plot strand bias with poisson test results
  plot = ggplot(strand_bias, aes( x = type,
                                  y = log2((strand_bias[,4]+0.1) /
                                             (strand_bias[,5]+0.1)),
                                  fill = type)) +
    scale_fill_manual(values = colors) +
    geom_bar(colour = "black", stat ="identity", position = "identity") +
    geom_text(
      aes(x = type,
          y = log2((strand_bias[,4]) / (strand_bias[,5]+0.1)),
          #ymax = log2((strand_bias[,4]) / (strand_bias[,5]+0.1)),
          label = significant,
          vjust = ifelse(sign(log2((strand_bias[,4]) /
                                     (strand_bias[,5]+0.1))) > 0, 0.5, 1)),
      size = 8,
      position = ggplot2::position_dodge(width = 1)) +
    theme_bw()  +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank()) +
    xlab("") +
    ylab(paste("log2(", var_names[1], "/", var_names[2], ")", sep = "") ) +
    scale_x_discrete(breaks=NULL)
  
  # If max_yaxis is TRUE, set max yaxis to maximum of all mutation types
  if (max_yaxis) { plot = plot + scale_y_continuous(limits = c(-max, max)) }
  
  # If different graphs for each mutation types is wanted, use method = "split"
  if (method == "split"){ plot = plot + facet_wrap(mutation ~ group, 
                                                   scales = "free_x",
                                                   nrow = length(levels(strand_bias$mutation)) )}
  else if (method == "combine"){ plot = plot + facet_wrap( ~ group, nrow = 1 )}
  return(plot)
}
