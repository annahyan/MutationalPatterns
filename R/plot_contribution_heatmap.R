#' Plot signature contribution heatmap
#' 
#' Plot relative contribution of signatures in a heatmap
#' 
#' @param contribution Signature contribution matrix
#' @param sig_order (Optional) Character vector with the desired order of the signature names for plotting. Optional.
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param cluster_samples (Optional) Hierarchically cluster samples based on eucledian distance. Default = TRUE.
#' @param cluster_mut_type (Optional) Signatures subset parameter for clustering
#' @param method (Optional) The agglomeration method to be used for hierarchical clustering. This should be one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
#' or "centroid" (= UPGMC).\cr
#' Default = "complete".
#' @param plot_values (Optional) Plot relative contribution values in heatmap.\cr
#' Default = F.
#' 
#' @return Heatmap with relative contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' 
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Set signature names as row names in the contribution matrix
#' rownames(nmf_res$contribution) = c("Signature A", "Signature B")
#' 
#' ## Define signature order for plotting
#' sig_order = c("Signature B", "Signature A")
#'
#'
#' ## Contribution heatmap with signatures in defined order
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"))
#' 
#' ## Contribution heatmap without sample clustering
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"), 
#'                           cluster_samples = FALSE, method = "complete")
#'
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}},
#' \code{\link{plot_contribution}}
#'
#' @export

# plotting function for relative contribution of signatures in heatmap
plot_contribution_heatmap = function(contribution, 
                                     sig_order, 
                                     type, 
                                     cluster_samples = TRUE, 
                                     cluster_mut_type, 
                                     method = "complete", 
                                     plot_values = FALSE)
{
  # check contribution argument
  if(class(contribution) == "list")
  { 
    for (m in names(contribution))
    {
        if (ncol(contribution[[m]]) == 1)
          cluster_samples = FALSE
    }
    combined = FALSE
  }
  else if (class(contribution) == "matrix")
  {
    if (ncol(contribution) == 1) cluster_samples = FALSE
    combined = TRUE
  } else {stop("contribution must be a named list")}
  
  if (!combined)
  {
    # Check type argument
    type = check_mutation_type(type)
    
    if (all(type %in% names(contribution))) 
        contribution = contribution[type]
    else {stop(paste("One or more values of 'type' is not found",
                     "in 'contribution'.",
                     "Run function without 'type' argument or give", 
                     "values which are in contribution list"))}
    
    # Hierarchically cluster samples
    if (cluster_samples)
    { 
      # Cluster signatures on mutation type
      if (missing(cluster_mut_type)) {cluster_mut_type = names(contribution)}
      else if (length(cluster_mut_type ) == 1)
      {
        if (cluster_mut_type == "all") 
          cluster_mut_type = c("snv","dbs","indel")
      }
      if (!all(cluster_mut_type %in% names(contribution)))
      {stop(paste("One or more values of 'cluster_mut_type' is not found", 
                  "in 'contribution'.", 
                  "Run function without 'cluster_mut_type' argument or give", 
                  "values which are in contribution list"))}
      
      cluster_mutations = c()
      
      # For each mutation type, get the mutation names for clustering
      for (m in cluster_mut_type)
      {
        if (any(grepl(m, names(contribution)))) 
          cluster_mutations = c(cluster_mutations, rownames(contribution[[m]]))
        else { cluster_mutations = cluster_mutations }
      }
      
      if (isEmpty(cluster_mutations))
        warning(paste("Values of 'cluster_mut_type' are not found in names of contribution list or do not match 'type'.",
                      "No clustering on mutation type is done"))
    } else {
      cluster_mutations = c()
    }
    
    contribution = do.call(rbind, contribution)
  } else 
  {
    cluster_mutations = c()
  }
  
  # check if there are signatures names in the contribution matrix
  if(is.null(row.names(contribution)))
    {stop("contribution must have row.names (signature names)")}
  # if no signature order is provided, check if signatures are clustered
  # if clustered, then clustered mutation types at front
  # else take same order as in the input matrix
  if(missing(sig_order))
  {
    if (!isEmpty(cluster_mutations))
    {
      sig_order = c(cluster_mutations, " ", 
                    rownames(contribution)[which(!(rownames(contribution) 
                                                   %in% cluster_mutations))])
      if (sig_order[length(sig_order)] == " ")
      {
        sig_order = sig_order[1:(length(sig_order)-1)] 
      }
    } else
      sig_order = rownames(contribution)
  } else {
    # check sig_order argument
    if(class(sig_order) != "character")
      {stop("sig_order must be a character vector")}
    if(length(sig_order) != nrow(contribution))
    {stop(paste("sig_order must have the same length as the number of",
                  "signatures in the contribution matrix"))}
    if(any(is.na(match(sig_order, row.names(contribution)))))
    {stop("sig_order must have the same signature names as in contribution")}
  }
    
  # transpose
  contribution = t(contribution)
  # relative contribution
  contribution_norm = contribution / rowSums(contribution)
  contribution_norm[which(contribution_norm == "NaN")] = 0
  
  # if cluster samples is TRUE, perform clustering
  if (cluster_samples)
  {
    # hiearchically cluster samples based on eucledian distance between 
    # relative contribution
    if (combined) {hc.sample = hclust(dist(contribution_norm), method = method)}
    else {hc.sample = hclust(dist(contribution_norm[,match(cluster_mutations, 
                                                           colnames(contribution_norm)),
                                                    drop=FALSE]),
                             method = method)}
    # order of samples according to hierarchical clustering
    sample_order = rownames(contribution)[hc.sample$order]
  } 
  else
  {
    sample_order = rownames(contribution)
  }

  Signature = NULL
  Sample = NULL
  Contribution = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL

  # melt data frame
  contribution_norm.m = melt(contribution_norm)
  # assign variable names
  colnames(contribution_norm.m) = c("Sample", "Signature", "Contribution")
  # change factor levels to the order for plotting
  contribution_norm.m$Sample = factor(contribution_norm.m$Sample, levels = sample_order)
  contribution_norm.m$Signature = factor(contribution_norm.m$Signature, levels = sig_order)
  
  if (cluster_samples)
  {
    if(!isEmpty(cluster_mutations) & length(cluster_mutations) != length(colnames(contribution)))
      contribution_norm.m = rbind(contribution_norm.m,
                                  data.frame("Sample" = unique(contribution_norm.m$Sample),
                                             "Signature" = " ",
                                             "Contribution" = NA))
  }
  
  # plot heatmap
  heatmap = ggplot(contribution_norm.m, aes(x=Signature, y=Sample, fill=Contribution, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \ncontribution", limits = c(0,1),
                         na.value = "white") +
    scale_x_discrete(breaks=sig_order, labels=sig_order) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x=NULL, y=NULL) 
  # if plot_values is TRUE, add values to heatmap
  if (plot_values)
  {
    heatmap = heatmap + geom_text(aes(label = round(Contribution, 2)), size = 3)
  }
  
  # if cluster_samples is TRUE, make dendrogram
  if (cluster_samples)
  {
    # get dendrogram
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse()+#expand = c(0.2, 0)) + 
      theme_dendro()
    # combine plots
    plot_final = cowplot::plot_grid(dendrogram, heatmap, align = 'h', rel_widths = c(0.3, 1), axis = "tb" )
  } 
  else
  {
    plot_final = heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(contribution_norm.m$Sample))))
  }

  return(plot_final)
}
