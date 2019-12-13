#' Plot strand per base substitution type
#' 
#' For each base substitution type and transcriptional strand the total number
#' of mutations and the relative contribution within a group is returned.
#' @param strand_counts A data.frame, result from strand_occurrences function
#' @param mode (Optional) Either "absolute" for absolute number of mutations, 
#' "relative" for relative contribution or "both" for both.\cr
#' Default = "relative"
#' @param colors (Optional) Named list with 6 value color vector "snv" for snv, 
#' 10 value color vector "dbs" for dbs.
#' For indels give same number of colors as there are classes.
#' @return Barplot
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
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

plot_strand = function(strand_counts, mode = "relative", colors)
{
  plots = list()
  
  if (class(strand_counts) == "list")
  { 
    strand_counts = do.call(rbind, strand_counts) 
    method = "split"
  } else { method = "combine" }
  
  # Set default colors
  if (missing(colors))
  {
    colors = c()
    if (any(grepl("snv", strand_counts$mutation))) { colors = c(colors, COLORS6) }
    if (any(grepl("dbs", strand_counts$mutation))) { colors = c(colors, COLORS10) }
    if (any(grepl("indel", strand_counts$mutation))) 
        colors = c(colors, COLORS_INDEL)
  }
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  type = NULL
  relative_contribution = NULL
  no_mutations = NULL
  
  if (mode == "both")
  {
    mode = "relative"
    mode_next = "absolute"
  } else {mode_next = "none"}
  
  strand_counts$mutation = factor(strand_counts$mutation, levels = unique(strand_counts$mutation))
  strand_counts$type = factor(strand_counts$type, levels = unique(strand_counts$type))
  
  # Plot relative contribution within each group
  if(mode == "relative")
  {
    plot = ggplot(strand_counts, aes(x=type,
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
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
    
    plots = c(plots, list(plot))
  }
  
  # Plot absolute contribution within each group
  if (mode == "absolute" | mode_next == "absolute")
  {
    plot = ggplot(strand_counts, aes(x=type,
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
      theme_bw() +
      scale_x_discrete(breaks=NULL) +
      xlab("")
    
      plots = c(plots, list(plot))
  }
  
  for (i in length(plots))
  {
    # If different graphs for each mutation types is wanted, use method = "split"
    if (method == "split"){ plots[[i]] = plots[[i]] + facet_wrap(mutation ~ group,  
                                                                 nrow = length(levels(strand_counts$mutation)),
                                                                 scales = "free_x")}
    else if (method == "combine"){ plots[[i]] = plots[[i]] + facet_wrap( ~ group,  nrow = 1 )}
  }
  
  plot = plot_grid(plotlist=plots)
  
  return(plot)
}