#' Plot signature contribution barplot
#' 
#' Plot contribution of signatures
#' 
#' @param contribution List of signature contribution matrices
#' @param signatures List of signature matrices
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param index (Optional) Sample subset parameter
#' @param coord_flip (Optional) Flip X and Y coordinates, default = FALSE
#' @param mode (Optional) Character stating to "relative" or "absolute" contribution of
#' mutations. Also possible to plot "both".\cr 
#' Default = "relative"
#' @param method (Optional) Character stating how to use the data. 
#' \itemize{
#'   \item{"split":} { Each mutation type has seperate count matrix}
#'   \item{"combine":} { Combined count matrix of all mutation types}
#' }   
#' Default is "split"
#' @param palette (Optional) A list of color palette like c("#FF0000", "#00FF00", "9999CC")
#' for each mutation type that will be used as colors in the plot.\cr
#' By default, ggplot2's colors are used to generate a palette.
#'
#' @return Stacked barplot with contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom reshape2 melt
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
                                type,
                                index=c(),
                                coord_flip=FALSE,
                                mode="relative",
                                method = "split",
                                palette=c())
{
    # check mode parameter
    if(!(mode == "relative" | mode == "absolute" | mode == "both"))
        stop("mode parameter should be either 'relative', 'absolute' or 'both'")
  
    # check mutation type
    if (missing(type)) type_default = TRUE
    else type_default = FALSE
    type = check_mutation_type(type)
    
    if (class(contribution) == "list")
    {
      if (!type_default)
      {
        if (all(type %in% names(contribution))) {contribution = contribution[type]}
        else {stop("One or more values of 'type' is not found in 'contribution'")}
      } else {
        type = names(contribution)
      }
      
      for (m in type)
      {
        if (is.null(rownames(contribution[[m]])))
          stop(paste("Provide contribution matrix for mutation type", 
                     m,
                     "with rownames for signatures"))
      }
      
      # optional subsetting if index parameter is provided
      if(length(index > 0)){
        for (m in type)
        {
          contribution[[m]] = contribution[[m]][,index]
        }
      }
      
      if (class(signatures) == "matrix"){
        for (m in names(contribution)){
          if (any(colnames(signatures) == rownames(contribution[[m]]))){
            signatures = list(signatures)
            names(signatures) = m
          }
        }
      }
    } else 
    {
      method = "combine"
      
      if (is.null(rownames(contribution)))
        stop("Provide contribution matrix with rownames for signatures")
      
      if (is.list(signatures))
      {
        keep_sigs = NULL
        for (m in names(signatures)){
          if (any(rownames(contribution) %in% colnames(signatures[[m]]))){
            keep_sigs <- c(keep_sigs,m)
          }
        }
        signatures = signatures[keep_sigs]
        signatures = do.call(rbind, signatures)
      }
    }
  
    if (length(palette) == 0)
    {
      if (class(signatures) == "matrix")
      {
        palette = default_colors_ggplot(ncol(signatures))
        names(palette) = colnames(signatures)  
      } else 
      {
        palette = list()
        for (m in names(signatures))
        {
          palette[[m]] = default_colors_ggplot(ncol(signatures[[m]]))
          names(palette[[m]]) = colnames(signatures[[m]])
        }
      }
    }
      
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    Sample = NULL
    Contribution = NULL
    Signature = NULL
    
    # Each mutation type has its own figure
    if (method == "split")
    {
      if (mode == "both")
      {
        mode = "relative"
        mode_next = "absolute"
      } else {mode_next = "none"}
      
      plots = list()
      
      for (m in type)
      {
        plots[[m]] = list()
        
        # Take all signatures with contribution more than 0
        sigs <- names(which(rowSums(contribution[[m]]) > 0))
        contribution[[m]] = contribution[[m]][which(rowSums(contribution[[m]]) > 0),
                                              ,drop=FALSE]
        
        if(length(sigs) == 1)
        {
          contribution[[m]] = t(as.matrix(contribution[[m]]))
          rownames(contribution[[m]]) = sigs
        }
        
        # Test if contribution is already relative
        if (all(round(colSums(contribution[[m]])) == 1))
        {
          warning(paste("Signature contributions are relative.",
                        "Plot absolute contribution is not possible"))
          mode = "relative"
          mode_next = "none"
        }
        
        if (mode == "relative")
        {
          # Plot contribution
          m_contribution = melt(contribution[[m]])
          colnames(m_contribution) = c("Signature", "Sample", "Contribution")
          m_contribution$Signature = factor(m_contribution$Signature, 
                                            levels = unique(m_contribution$Signature))
          
          plot = ggplot(m_contribution,
                        aes(x = factor(Sample),
                            y = Contribution,
                            fill = Signature,
                            order = Sample)) +
            geom_bar(position = "fill", stat="identity")  +
            # ylabel
            labs(x = "", y = "Relative contribution") +
            # white background
            theme_bw() +
            # default or custom palette
            scale_fill_manual(name="Signature", values = palette[[m]]) +
            # no gridlines
            theme(panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank()) +
            theme(panel.grid.minor.y=element_blank(),
                  panel.grid.major.y=element_blank())
          
          plots[[m]][["relative"]] = plot
        }
      
        # Handle the absolute mode.
        if (mode == "absolute" | mode_next == "absolute")
        {
          if(missing(signatures))
            stop(paste("For contribution plotting in mode 'absolute':",
                       "also provide signatures matrix"))

          # total number of mutations per signature
          total_signatures = colSums(signatures[[m]])[which(colnames(signatures[[m]]) %in% rownames(contribution[[m]]))] 
          
          # calculate signature contribution in absolute number of signatures
          abs_contribution = contribution[[m]] * total_signatures
          
          # Plot contribution
          m_contribution = melt(abs_contribution)
          colnames(m_contribution) = c("Signature", "Sample", "Contribution")
          m_contribution$Signature = factor(m_contribution$Signature, 
                                            levels = unique(m_contribution$Signature))
          
          plot = ggplot(m_contribution, aes(x = factor(Sample),
                                            y = Contribution,
                                            fill = factor(Signature),
                                            order = Sample)) + 
            geom_bar(stat="identity")  +  
            # ylabel
            labs(x = "", y = "Absolute contribution \n (no. mutations)") +  
            # white background
            theme_bw() +
            # default or custom palette
            scale_fill_manual(name="Signature", values = palette[[m]]) +
            # no gridlines
            theme(panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank()) +
            theme(panel.grid.minor.y=element_blank(),
                  panel.grid.major.y=element_blank())
          
          plots[[m]][["absolute"]] = plot
        }
      
        if (mode_next == "none")
        {
          # Handle coord_flip.
          if (coord_flip)
            plots[[m]][[1]] = plots[[m]][[1]] + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
          else
            plots[[m]][[1]] = plots[[m]][[1]] + xlim(levels(factor(m_contribution$Sample)))
          
          plots[[m]] = plots[[m]][[1]]

        } else 
        {
          for (n in 1:length(plots[[m]]))
          {
            # Handle coord_flip.
            if (coord_flip)
              plots[[m]][[n]] = plots[[m]][[n]] + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
            else
              plots[[m]][[n]] = plots[[m]][[n]] + xlim(levels(factor(m_contribution$Sample)))
          }
          
          plots[[m]][[1]] = plots[[m]][[1]] + theme(legend.title = element_blank())
        }
      }
      
      # Make list out of list of lists  
      plotlist = list()
      for (m in names(plots))
      {
        if (length(class(plots[[m]])) == 1)
        {
          plotlist = c(plotlist, list(plots[[m]][[1]]))
          plotlist = c(plotlist, list(plots[[m]][[2]]))
        } else {
          plotlist = c(plotlist, list(plots[[m]]))
        }
      }
      
      if (mode_next == "absolute") 
        return(plot_grid(plotlist = plotlist, ncol=2, align="hv", axis="tblr"))
      else 
        return(plot_grid(plotlist = plotlist, ncol=1, align="v", axis = "lr"))
    } else if (method == "combine")
    {
      if (mode == "both")
      {
        mode = "relative"
        mode_next = "absolute"
      } else {mode_next = "none"}
      
      if (class(contribution) == "list")
      {
        colsums = c()
        total_signatures = list()
        abs_contribution = list()
        
        for (m in type){
          if (all(round(colSums(contribution[[m]])) ==1 ))
            colsums = c(colsums, TRUE)
          else 
            colsums = c(colsums, FALSE)
          
          # total number of mutations per signature
          total_signatures[[m]] = colSums(signatures[[m]])[which(colnames(signatures[[m]]) %in% rownames(contribution[[m]]))] 
          
          # calculate signature contribution in absolute number of signatures
          abs_contribution[[m]] = contribution[[m]] * total_signatures[[m]]
        }
        
        # Test if contribution is already relative
        if (all(colsums == TRUE))
        {
          warning(paste("Signature contributions are relative.",
                        "Plot absolute contribution is not possible"))
          mode = "relative"
          mode_next = "none"
        }
        
        names(contribution) = NULL
        contribution = do.call(rbind, contribution)
        names(abs_contribution) = NULL
        abs_contribution = do.call(rbind, abs_contribution)
      } else {
        # total number of mutations per signature
        total_signatures = colSums(signatures)[which(colnames(signatures) %in% rownames(contribution))]  
        
        # calculate signature contribution in absolute number of signatures
        abs_contribution = contribution * total_signatures
      }
      
      # Take all signatures with contribution more than 0
      contribution = contribution[which(rowSums(contribution) > 0),
                                  ,drop=FALSE]
      abs_contribution = abs_contribution[which(rowSums(abs_contribution) > 0),
                                          ,drop=FALSE]
      
      palette = unlist(unname(palette))
      
      plots = list()
      
      if (mode == "relative")
      {
          # Plot contribution
          m_contribution = melt(contribution)
          colnames(m_contribution) = c("Signature", "Sample", "Contribution")
          m_contribution$Signature = factor(m_contribution$Signature, 
                                            levels = unique(m_contribution$Signature))
  
          plot = ggplot(m_contribution,
                          aes(x = factor(Sample),
                              y = Contribution,
                              fill = Signature,
                              order = Sample)) +
              geom_bar(position = "fill", stat="identity")  +
              # ylabel
              labs(x = "", y = "Relative contribution") +
              # white background
              theme_bw() +
              # default or custom palette
              scale_fill_manual(name="Signature", values = palette) +
              # no gridlines
              theme(panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank()) +
              theme(panel.grid.minor.y=element_blank(),
                      panel.grid.major.y=element_blank())
          
          plots = c(plots, list(plot))
      }
  
      # Handle the absolute mode.
      if (mode == "absolute" | mode_next == "absolute")
      {
          if(missing(signatures))
              stop(paste("For contribution plotting in mode 'absolute':",
                          "also provide signatures matrix"))
          
          # Plot contribution
          m_contribution = melt(abs_contribution)
          colnames(m_contribution) = c("Signature", "Sample", "Contribution")
          m_contribution$Signature = factor(m_contribution$Signature, 
                                            levels = unique(m_contribution$Signature))
  
          plot = ggplot(m_contribution, aes(x = factor(Sample),
                                              y = Contribution,
                                              fill = Signature,
                                              order = Sample)) + 
              geom_bar(stat="identity")  +  
              # ylabel
              labs(x = "", y = "Absolute contribution \n (no. mutations)") +  
              # white background
              theme_bw() +
              # default or custom palette
              scale_fill_manual(name="Signature", values = palette) +
              # no gridlines
              theme(panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank()) +
              theme(panel.grid.minor.y=element_blank(),
                      panel.grid.major.y=element_blank())
          
          plots = c(plots, list(plot))
      }
      
      if (mode_next == "none")
      {
        # Handle coord_flip.
        if (coord_flip)
            plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
        else
            plot = plot + xlim(levels(factor(m_contribution$Sample)))
      } else 
      {
        for (m in 1:length(plots))
        {
          # Handle coord_flip.
          if (coord_flip)
            plots[[m]] = plots[[m]] + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
          else
            plots[[m]] = plots[[m]] + xlim(levels(factor(m_contribution$Sample)))
        }
        
        plot = plot_grid(plotlist=plots,ncol=1, align="v")
      }
      
      return(plot)
  }
}
