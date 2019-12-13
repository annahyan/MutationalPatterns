#' Compare two mutational profiles
#'
#' Plots two mutational profiles and their difference, reports the residual
#' sum of squares (RSS).
#'
#' @param profile1 First mutational profile
#' @param profile2 Second mutational profile
#' @param profile_names (Optional) Character vector with names of the mutations profiles
#' used for plotting.\cr
#' Default = c("profile 1", "profile 2")
#' @param profile_ymax (Optional) Maximum value of y-axis (relative contribution) for
#' profile plotting.\cr
#' Default = 0.2
#' @param diff_ylim (Optional) Y-axis limits for profile difference plot.\cr
#' Default = c(-0.02, 0.02)
#' @param colors (Optional) List of color vectors with same length as mutational classes per type:\cr
#' 6 for SNV, 10 for DBS and number of indel classes for wanted context (6 for "predefined" and 16 for "cosmic").\cr
#' Default colors are predefined in MutationalPatterns
#' @param condensed (Optional) More condensed plotting format. Default = F.
#' @return mutational profile plot of profile 1, profile 2 and their difference
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
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
#' ## Compare the reconstructed 96-profile of sample 1 with orignal profile
#' plot_compare_profiles(mut_mat[,1],
#'                         nmf_res$reconstructed[,1],
#'                         profile_names = c("Original", "Reconstructed"))
#'
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{extract_signatures}}
#'
#' @export


plot_compare_profiles = function(profile1,
                                    profile2,
                                    profile_names = c("profile 1", "profile 2"),
                                    profile_ymax = 0.2,
                                    diff_ylim = c(-0.02, 0.02),
                                    colors,
                                    condensed = FALSE)
{
    # Check if both profiles are from the same mutation type
    if (any(is.na(match(names(profile1), names(profile2)))) | any(is.na(match(names(profile1), names(profile2)))))
    { stop("Mutations of profiles do not match. Is the same mutation type given?")}
    
    # Find the mutation type(s) for coloring and plotting the contexts
    if (all(names(profile1) %in% TRIPLETS_96)) { type = "snv" }
    else if (all(names(profile1) %in% DBS)) { type = "dbs" }
    else if (all(names(profile1) %in% c(TRIPLETS_96, DBS))) { type = c("snv", "dbs")}
    else 
    {
        if (all(names(profile1) %in% INDEL_CONTEXT)) { type = "indel" }
        else if (all(names(profile1) %in% c(TRIPLETS_96, INDEL_CONTEXT)))
        {
            warning("Mutation type of profile1 is unknown. Treated as combined mutation type")
            type = c("snv", "indel")
        }
        else if (all(names(profile1) %in% c(DBS, INDEL_CONTEXT)))
        {
            warning("Mutation type of profile1 is unknown. Treated as combined mutation type")
            type = c("dbs", "indel")
        }
        else if (all(names(profile1) %in% c(TRIPLETS_96, DBS, INDEL_CONTEXT)))
        {
            warning("Mutation type of profile1 is unknown. Treated as combined mutation type")
            type = c("snv", "dbs", "indel")
        } else {
            stop("Mutations in profile 1 are not found in preset SNV and DBS or in given INDEL context")
        }
    }
  
    # if colors parameter not provided, set to default colors
    if(missing(colors))
    {
      colors = list()
      if ("snv" %in% type) { colors[["snv"]] = COLORS6 }
      if ("dbs" %in% type) { colors[["dbs"]] = COLORS10 }
      if ("indel" %in% type) { colors[["indel"]] = COLORS_INDEL }
    }
  
    # Get relative profiles and difference
    s1_relative = profile1 / sum(profile1)
    s2_relative = profile2 / sum(profile2)
    diff = s1_relative - s2_relative

    # residual sum of squares
    RSS = sum(diff^2)
    RSS = format(RSS, scientific = TRUE, digits = 3)
    
    # calculate cosine similarity between the two profiles
    cosine_sim = cos_sim(profile1, profile2)
    # round
    cosine_sim = round(cosine_sim, 3)
    
    x = cbind(s1_relative, s2_relative, diff)
    colnames(x) = c(profile_names, "Difference")
    
    # Get context and substitutions info of mutation types
    substitutions = list()
    context = list()
    
    if ("snv" %in% type)
    {
      substitutions[["snv"]] = SUBSTITUTIONS_96
      index = c(rep(1,1,16), rep(2,1,16), rep(3,1,16),
                rep(4,1,16), rep(5,1,16), rep(6,1,16))
      # Context
      context_snv = CONTEXTS_96
      
      # Replace mutated base with dot
      substring(context_snv,2,2) = "."
      context[["snv"]] = context_snv
    } 
    if ("dbs" %in% type)
    {
      substitutions[["dbs"]] = unlist(lapply(as.list(SUBSTITUTIONS_DBS), function(sub)
      {
        sub <- unlist(strsplit(sub, ">"))[1]
        l <- length(which(startsWith(DBS, sub)))
        return(rep(paste0(sub,">NN"), l))
      }))
      context[["dbs"]] = ALT_DBS
    }
    if ("indel" %in% type)
    {
      substitutions[["indel"]] = INDEL_CLASS
      context[["indel"]] = do.call(rbind, strsplit(INDEL_CONTEXT, "\\."))[,lengths(strsplit(INDEL_CONTEXT, "\\."))[1]]
    }
    
    # Translate lists into vector
    colors = unname(unlist(colors))
    substitutions = unname(unlist(substitutions))
    context = unname(unlist(context))
    
    # Construct dataframe for plotting
    df = data.frame(substitution = substitutions, context = context)
    rownames(x) = NULL
    df2 = cbind(df, as.data.frame(x))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    df3$substitution <- factor(df3$substitution, levels = unique(df3$substitution))

    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    value = NULL
    substitution = NULL
    Sample = NULL
    Contribution = NULL
    Signature = NULL

    # Add dummy non_visible data points to force y axis limits per facet
    df4 = data.frame(substitution = rep(substitutions[1], 4),
                        context = rep(context[1],4),
                        variable = c(profile_names, "Difference", "Difference"),
                        value = c(profile_ymax,
                                    profile_ymax,
                                    diff_ylim[1],
                                    diff_ylim[2]))
    if (condensed)
    {
      plot = ggplot(data=df3, aes(x=context,
                                  y=value,
                                  fill=substitution,
                                  width=1)) +
        geom_bar(stat="identity",
                 position = "identity",
                 colour="black", size=.2) +
        geom_point(data = df4, aes(x = context,
                                   y = value), alpha = 0) +
        scale_fill_manual(values=colors) +
        facet_grid(variable ~ substitution, scales = "free") +
        ylab("Relative contribution") +
        # ylim(-yrange, yrange) +
        # no legend
        guides(fill=FALSE) +
        # white background
        theme_bw() +
        ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = "")) +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
              axis.text.y=element_text(size=8),
              axis.title.x=element_text(size=12),
              axis.text.x=element_text(size=5,angle=90,vjust=0.4),
              strip.text.x=element_text(size=14),
              strip.text.y=element_text(size=14),
              panel.grid.major.x = element_blank(),
              panel.spacing.x = unit(0, "lines"))
      
    } else {
    plot = ggplot(data=df3, aes(x=context,
                                y=value,
                                fill=substitution,
                                width=0.6)) +
        geom_bar(stat="identity",
                    position = "identity",
                    colour="black", size=.2) +
        geom_point(data = df4, aes(x = context,
                                    y = value), alpha = 0) +
        scale_fill_manual(values=colors) +
        facet_grid(variable ~ substitution, scales = "free_y") +
        ylab("Relative contribution") +
        # ylim(-yrange, yrange) +
        # no legend
        guides(fill=FALSE) +
        # white background
        theme_bw() +
        ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = "")) +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=12),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=14),
                strip.text.y=element_text(size=14),
                panel.grid.major.x = element_blank())
    }
    return(plot)
}
