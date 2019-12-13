#' Plot signature strand bias
#' 
#' Plot strand bias per mutation type for SNV, DBS and indels.
#' 
#' @param signatures_strand_bias Named list of signatures matrices
#' @param colors (Optional) Named list with 6 value color vector "snv" for snv, 10 value color vector "dbs" for dbs.
#' For indels give same number of colors as there are classes. \cr
#' 'native' contains 6 classes and 'cosmic' contains 16 classes
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @return Barplot
#'
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom plyr adply
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
#'
#' nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
#'                                         package="MutationalPatterns"))
#'
#' ## Provide column names for the plot.
#' colnames(nmf_res_strand$signatures) = c("Signature A", "Signature B")
#'
#' plot_signature_strand_bias(nmf_res_strand$signatures) 
#'
#' @seealso
#' \code{link{extract_signatures}},
#' \code{link{mut_matrix()}}
#'
#' @export

plot_signature_strand_bias = function(signatures_strand_bias, type, colors)
{
    # Check mutation type argument
    if (missing(type)) type_default = TRUE
    else type_default = FALSE
    type = check_mutation_type(type)  
    
    if (class(signatures_strand_bias) == "list")
    {  
      method = "split"
    } else { method = "combine" }
    
    # If strand bias is a matrix, find the mutation types present by 
    # comparing the rownames with defined contexts
    if (class(signatures_strand_bias) == "matrix")
    {
      types_info = do.call(rbind, strsplit(rownames(signatures_strand_bias), "-"))
      types = types_info[,1]
      strand = unique(types_info[,2])
      
      if (all(types %in% TRIPLETS_96)){
        signatures_strand_bias = list("snv"=signatures_strand_bias)
      } else if (all(types %in% DBS))
      {
        signatures_strand_bias = list("dbs"=signatures_strand_bias)
      } else if (all(types %in% c(TRIPLETS_96,DBS)))
      {
        signatures_strand_bias = list("snv"=signatures_strand_bias[types %in% TRIPLETS_96,],
                                      "dbs"=signatures_strand_bias[types %in% DBS,])
        method = "combine"
      } else if (all(types %in% INDEL_CONTEXT))
      {
        signatures_strand_bias = list("indel"=signatures_strand_bias)
      } else if (all(types %in% c(TRIPLETS_96,INDEL_CONTEXT)))
      {
        signatures_strand_bias = list("snv"=signatures_strand_bias[types %in% TRIPLETS_96,],
                          "indel"=signatures_strand_bias[types %in% INDEL_CONTEXT,])
        method = "combine"
      } else if (all(types %in% c(DBS, INDEL_CONTEXT)))
      {
        signatures_strand_bias = list("dbs"=signatures_strand_bias[types %in% DBS,],
                          "indel"=signatures_strand_bias[types %in% INDEL_CONTEXT,])
        method = "combine"
      } else if (all(types %in% c(TRIPLETS_96, DBS, INDEL_CONTEXT)))
      {
        signatures_strand_bias = list("snv"=signatures_strand_bias[types %in% TRIPLETS_96,],
                          "dbs"=signatures_strand_bias[types %in% DBS,],
                          "indel"=signatures_strand_bias[types %in% INDEL_CONTEXT,])
        method = "combine"
      } else {
        stop("Mutation matrix is not a list and mutation types could not be derived")
      }
    }
    
    if (type_default) type = names(signatures_strand_bias)
    else type = intersect(type, names(signatures_strand_bias))
    
    if (missing(colors))
    {
      colors = list()
      if ("snv" %in% names(signatures_strand_bias)) { colors[["snv"]] = COLORS6 }
      if ("dbs" %in% names(signatures_strand_bias)) { colors[["dbs"]] = COLORS10 }
      if ("indel" %in% names(signatures_strand_bias)) { colors[["indel"]] = COLORS_INDEL }
    }
  
    # Set strand info
    strand = list("snv" = rep(c("U", "T"), 96), 
                  "dbs" = rep(c("U", "T"), 78),
                  "indel" = rep(c("U", "T"), length(INDEL_CONTEXT)))
    
    # Set substitutions of mutation types
    substitutions = list("snv" = list("type"=SUBSTITUTIONS, 
                                      "strand"=rep(SUBSTITUTIONS, each = 32)),
                         "dbs" = list("type"=SUBSTITUTIONS_DBS, 
                                      "strand"=rep(do.call(rbind, strsplit(DBS, ">"))[,1], each = 2)),
                         "indel" = list("type"=unique(paste(INDEL_CLASS_HEADER,INDEL_CLASS, sep=".")),
                                        "strand"=rep(paste(INDEL_CLASS_HEADER,INDEL_CLASS, sep="."), each=2)))
    
    plots = list()
    
    if (method == "split"){
      
      # Get all info for the mutation types asked for
      strand = strand[type]
      substitutions = substitutions[type]
      signatures_strand_bias = signatures_strand_bias[type]
      
      # Test if strand bias is for replication of transcription
      if (endsWith(rownames(signatures_strand_bias[[1]])[1], "left") | endsWith(rownames(signatures_strand_bias[[1]])[1], "right"))
      {
        mode = "replication" 
      } else if (endsWith(rownames(signatures_strand_bias[[1]])[1], "transcribed") | endsWith(rownames(signatures_strand_bias[[1]])[1], "untranscribed"))
      { 
        mode = "transcription"
      }
      
      for (m in names(signatures_strand_bias))
      {
        # aggregate by strand and type
        sum_per_type = stats::aggregate(signatures_strand_bias[[m]],
                                        by=list(strand[[m]], substitutions[[m]][["strand"]]),
                                        FUN=sum)
        
        sum_per_strand = stats::aggregate(signatures_strand_bias[[m]],
                                          by=list(strand[[m]]),
                                          FUN=sum)
        
        # melt data frames
        sum_per_strand =  melt(sum_per_strand)
        colnames(sum_per_strand) = c("strand", "Signature", "value")
        sum_per_type =  melt(sum_per_type)
        colnames(sum_per_type) = c("strand", "type", "Signature", "value")
        
        # ratio per signature per type
        ratio = as.matrix(subset(sum_per_type, strand == "T")$value /
                            subset(sum_per_type, strand == "U")$value)
        
        ratio_per_type_per_signature = cbind(subset(sum_per_type,
                                                    strand == "T")[,2:3],
                                             ratio)
        
        # binomial test per type per signature
        size = c()
        observed = c()
        transcribed = c()
        untranscribed = c()
        for (s in unique(sum_per_type$Signature))
        {
          for (t in unique(sum_per_type$type))
          {
            sub = subset(sum_per_type, Signature==s & type==t)
            size = c(size, sum(sub$value))
            observed = c(observed, subset(sub, strand == "T")$value)
            transcribed = c(transcribed, subset(sub, strand == "T")$value)
            untranscribed = c(untranscribed, subset(sub, strand == "U")$value)
          }
        }
        
        # names of signatures
        signatures = colnames(signatures_strand_bias[[m]])
        
        # No. signatures
        n = length(signatures)
        
        # Observed: observed no. mutations on transcribed strand rounded to integer
        # Combine counts in one data.frame
        if (m == "snv")
        {  
          stats_per_type = data.frame(
            Signature = rep(signatures, each=6),
            type = rep(substitutions[[m]][["type"]],n),
            size = as.integer(size),
            transcribed = transcribed,
            untranscribed = untranscribed,
            observed = as.integer(observed))
        } else if (m == "dbs")
        {
          stats_per_type = data.frame(
            Signature = rep(signatures, each=10),
            type = rep(substitutions[[m]][["type"]],n),
            size = as.integer(size),
            transcribed = transcribed,
            untranscribed = untranscribed,
            observed = as.integer(observed))
        } else if (m == "indel")
        {
          stats_per_type = data.frame(
            Signature = rep(signatures, each=16),
            type = rep(substitutions[[m]][["type"]],n),
            size = as.integer(size),
            transcribed = transcribed,
            untranscribed = untranscribed,
            observed = as.integer(observed))
        }
        # Perform binomial test
        stats_per_type = plyr::adply(
          stats_per_type,
          1,
          function(x) binomial_test(0.5, x$size, x$observed))
        
        # Calculate ratio 
        ratio_per_type_per_signature = cbind(ratio_per_type_per_signature,
                                             stats_per_type)
        
        strand_bias_per_type_df = melt(ratio_per_type_per_signature[,c(1,2,3,12)])
        
        # Find maximum y value for plotting, round up to integer
        maximum = ceiling(max(abs(log2(strand_bias_per_type_df$value))))
        
        # These variables will be available at run-time, but not at compile-time.
        # To avoid compiling trouble, we initialize them to NULL.
        Signature = NULL
        type = NULL
        value = NULL
        significant = NULL
        
        # Plot
        plot = ggplot(strand_bias_per_type_df,
                      aes(x=type, y=log2(value), fill=type)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          scale_fill_manual(values=colors[[m]]) +
          scale_y_continuous(limits = c(-maximum, maximum)) +
          facet_grid(Signature ~ .) +
          theme_bw() + 
          scale_x_discrete(breaks=NULL) +
          xlab("") +
          geom_text(
            aes(x = type,
                y = log2(value),
                ymax = log2(value), 
                label = significant,
                vjust = ifelse(sign(log2(value)) > 0, 0.5, 1)), 
            size = 8, position = ggplot2::position_dodge(width=1))
        
        if (mode == "replication") { plot = plot + ylab("log2(left/right)")}
        else if (mode == "transcription") { plot = plot + ylab("log2(transcribed/untranscribed)")}
        
        plots[[m]] = plot
      }
      
      return(plot_grid(plotlist = plots, ncol = 1, align = "v"))
    } else if (method == "combine")
    {
      # Get all info for the mutation types asked for
      strand = strand[type]
      substitutions = substitutions[type]
      signatures_strand_bias = signatures_strand_bias[type]
      
      # Unlist all data
      strand = unlist(unname(strand))
      substitutions_new = list("type"=c(), 
                               "strand"=c())
      for (m in names(substitutions))
      {
        substitutions_new[["type"]] = c(substitutions_new$type, 
                                        substitutions[[m]]$type)        
        substitutions_new[["strand"]] = c(substitutions_new$strand,
                                          substitutions[[m]]$strand)
      }
      substitutions = substitutions_new
      
      signatures_strand_bias = do.call(rbind, signatures_strand_bias)
      
      # Test if strand bias is for replication or transcription
      if (endsWith(rownames(signatures_strand_bias)[1], "left") | endsWith(rownames(signatures_strand_bias)[1], "right"))
      {
        mode = "replication" 
      } else if (endsWith(rownames(signatures_strand_bias)[1], "transcribed") | endsWith(rownames(signatures_strand_bias)[1], "untranscribed"))
      { 
        mode = "transcription"
      }

      # aggregate by strand and type
      sum_per_type = stats::aggregate(signatures_strand_bias,
                                      by=list(strand, substitutions[["strand"]]),
                                      FUN=sum)
      
      sum_per_strand = stats::aggregate(signatures_strand_bias,
                                        by=list(strand),
                                        FUN=sum)
      
      # melt data frames
      sum_per_strand =  melt(sum_per_strand)
      colnames(sum_per_strand) = c("strand", "Signature", "value")
      sum_per_type =  melt(sum_per_type)
      colnames(sum_per_type) = c("strand", "type", "Signature", "value")
      
      # ratio per signature per type
      ratio = as.matrix(subset(sum_per_type, strand == "T")$value /
                          subset(sum_per_type, strand == "U")$value)
      
      ratio_per_type_per_signature = cbind(subset(sum_per_type,
                                                  strand == "T")[,2:3],
                                           ratio)
    
      # binomial test per type per signature
      size = c()
      observed = c()
      transcribed = c()
      untranscribed = c()
      for (s in unique(sum_per_type$Signature))
      {
        for (t in unique(sum_per_type$type))
        {
          sub = subset(sum_per_type, Signature==s & type==t)
          size = c(size, sum(sub$value))
          observed = c(observed, subset(sub, strand == "T")$value)
          transcribed = c(transcribed, subset(sub, strand == "T")$value)
          untranscribed = c(untranscribed, subset(sub, strand == "U")$value)
        }
      }
      
      # names of signatures
      signatures = colnames(signatures_strand_bias)
      
      # No. signatures
      n = length(signatures)
      
      # Observed: observed no. mutations on transcribed strand rounded to integer
      # Combine counts in one data.frame
      stats_per_type = data.frame(
        Signature = rep(signatures, each=length(substitutions[["type"]])),
        type = rep(substitutions[["type"]],n),
        size = as.integer(size),
        transcribed = transcribed,
        untranscribed = untranscribed,
        observed = as.integer(observed))
      
      # Perform binomial test
      stats_per_type = plyr::adply(
        stats_per_type,
        1,
        function(x) binomial_test(0.5, x$size, x$observed))
      
      # Calculate ratio 
      ratio_per_type_per_signature = cbind(ratio_per_type_per_signature,
                                           stats_per_type)
      
      strand_bias_per_type_df = melt(ratio_per_type_per_signature[,c(1,2,3,12)])
      strand_bias_per_type_df$type = factor(strand_bias_per_type_df$type, levels = unique(substitutions[["strand"]]))
      
      # Find maximum y value for plotting, round up to integer
      maximum = ceiling(max(abs(log2(strand_bias_per_type_df$value))))
      
      # These variables will be available at run-time, but not at compile-time.
      # To avoid compiling trouble, we initialize them to NULL.
      Signature = NULL
      type = NULL
      value = NULL
      significant = NULL
      
      # Plot
      plot = ggplot(strand_bias_per_type_df,
                    aes(x=type, y=log2(value), fill=type)) +
        geom_bar(stat="identity", position="dodge", color="black") +
        scale_fill_manual(values=unname(unlist(colors))) +
        scale_y_continuous(limits = c(-maximum, maximum)) +
        facet_grid(Signature ~ .) +
        theme_bw() + 
        scale_x_discrete(breaks=NULL) +
        xlab("") +
        geom_text(
          aes(x = type,
              y = log2(value),
              ymax = log2(value), 
              label = significant,
              vjust = ifelse(sign(log2(value)) > 0, 0.5, 1)), 
          size = 8, position = ggplot2::position_dodge(width=1))
      
      if (mode == "replication") { plot = plot + ylab("log2(left/right)") }
      else if (mode == "transcription") { plot = plot + ylab("log2(transcribed/untranscribed)")}
      
      return(plot)
    }
}
