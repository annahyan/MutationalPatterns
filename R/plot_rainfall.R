#' Plot genomic rainfall 
#'
#' Rainfall plot visualizes the types of mutations and intermutation distance
#' @details
#' Rainfall plots can be used to visualize the distribution of mutations
#' along the genome or a subset of chromosomes. The distance of a mutation
#' with the mutation prior to it (the intermutation distance) is plotted on
#' the y-axis on a log scale.
#'
#' The colour of the points indicates the base substitution type.
#' Clusters of mutations with lower intermutation distance represent mutation
#' hotspots.
#'
#' @param vcf CollapsedVCF object
#' @param chromosomes Vector of chromosome/contig names of the reference
#' genome to be plotted
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param method (Optional) Character stating how to plot the results. method = "split" will give 
#' seperate plots for each mutation type, whereas method = "combine" will give one plot with all
#' mutation types.\cr
#' Default is "split"
#' @param title (Optional) Plot title
#' @param colors (Optional) Named list with 6 value color vector "snv" for snv, 
#' 10 value color vector "dbs" for dbs.
#' For indels give same number of colors as there are classes. \cr
#' @param cex (Optional) Point size, default = 2.5
#' @param cex_text (Optional) Text size, default = 3
#' @param ylim (Optional) Maximum y value (genomic distance), default = 1e+8
#' @return Rainfall plot
#'
#' @import ggplot2
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqnames
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#' 
#' # Specify chromosomes of interest.
#' chromosomes = names(genome(vcfs[[1]])[1:22])
#'
#' ## Do a rainfall plot for all chromosomes:
#' plot_rainfall(vcfs[[1]],
#'                 title = names(vcfs[1]),
#'                 chromosomes = chromosomes,
#'                 cex = 1)
#'
#' ## Or for a single chromosome (chromosome 1):
#' plot_rainfall(vcfs[[1]],
#'                 title = names(vcfs[1]),
#'                 chromosomes = chromosomes[1],
#'                 cex = 2)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

plot_rainfall <- function(vcf, chromosomes, type, method = "split", title = "", colors, cex = 2.5,
                            cex_text = 3, ylim = 1e+08)
{
    # Check the mutation type argument
    type = check_mutation_type(type)
    
    # If colors parameter not provided, set to default colors
    if(missing(colors))
    { 
      colors_list=list("snv"=COLORS6, 
                       "dbs"=COLORS10,
                       "indel"=COLORS_INDEL)
    } else { colors_list = colors }
    
    # Check color vector length
    if ("snv" %in% type)
      if (length(colors_list$snv) != 6)
        stop("Provide colors vector for single base substitutions with length 6")
    if ("dbs" %in% type)
      if (length(colors_list$dbs) != 10)
        stop("Provide colors vector for double base substitutions with length 10")
    
    indel_color_number = 1
    for (i in 2:length(INDEL_CLASS))
    {
      if (INDEL_CLASS[i-1] != INDEL_CLASS[i]) { indel_color_number = indel_color_number + 1 }
    }
    
    if(length(unique(colors_list$indel)) != indel_color_number)
      stop("Provide indel colors vector with length same number of classes")

    # Get all the substitutions for each mutation type  
    substitutions_list = list("snv"=SUBSTITUTIONS, 
                              "dbs"=SUBSTITUTIONS_DBS, 
                              "indel"=unique(paste(INDEL_CLASS_HEADER,INDEL_CLASS, sep=".")))
    
    # For each mutation type, get all mutations from the vcf
    muts = list()
    for (m in type)
    {
      if (m == "snv") { muts[[m]] = which(nchar(as.character(vcf$REF))==1 & nchar(as.character(unlist(vcf$ALT)))==1) }
      else if (m == "dbs") { muts[[m]] = which(nchar(as.character(vcf$REF))==2 & nchar(as.character(unlist(vcf$ALT)))==2) }
      else if (m == "indel") { muts[[m]] = which(nchar(as.character(vcf$REF)) != nchar(as.character(unlist(vcf$ALT))) & 
                                                   (nchar(as.character(vcf$REF)) == 1 | nchar(as.character(unlist(vcf$ALT))) ==1 )) }
    }
    
    if (method == "split")
    {
      plots = list()
      
      for (mut in type)
      {
        input_vcf = vcf[which(1:length(vcf) %in% muts[[mut]]),]
        input_vcf = input_vcf[order(input_vcf),]
        
        colors = colors_list[[mut]]
      
        # get chromosome lengths of reference genome
        chr_length = seqlengths(input_vcf)
        
        # subset
        chr_length = chr_length[names(chr_length) %in% chromosomes]
        
        # cumulative sum of chromosome lengths
        chr_cum = c(0, cumsum(as.numeric(chr_length)))
        
        # Plot chromosome labels without "chr"
        names(chr_cum) = names(chr_length)
        labels = gsub("chr", "", names(chr_length))
        
        # position of chromosome labels
        m=c()
        for(i in 2:length(chr_cum))
          m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)
        
        # mutation characteristics
        types = loc = dist = chrom = c()
        
        # for each chromosome get types, location, distance and chromosome information
        # of all mutations
        for(i in 1:length(chromosomes))
        {
          chr_subset = input_vcf[seqnames(input_vcf) == chromosomes[i]]
          
          n = length(chr_subset)
          if(n<=1){next}
          if (mut == "snv")
            type_list = mut_type(chr_subset, type = mut)
          else if (mut == "dbs")
          {
            type_list = mut_type(chr_subset, type = mut)
            for (j in 1:length(type_list)) type_list[j] = paste0(substr(type_list[j],1,2),">NN")
          } else if (mut == "indel")
          {
            type_list = mut_context(chr_subset, ref_genome, type = mut)
            type_list = do.call(rbind, strsplit(type_list, "\\."))[,c(1,2,4)]
            type_list = list("indel"=paste(type_list[,1], type_list[,2], type_list[,3], sep="."))
          }
  
          types = c(types, unname(unlist(type_list))[-1])
          loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
          dist = c(dist, diff(start(chr_subset)))
          chrom = c(chrom, rep(chromosomes[i],n-1))
        }
        
        data = data.frame(type = types,
                          location = loc,
                          distance = dist,
                          chromosome = chrom)
        
        # Removes colors based on missing mutation types.  This prevents colors from
        # shifting when comparing samples with low mutation counts.
        substitutions = substitutions_list[[mut]]
        typesin = substitutions %in% data$type
        colors = colors[typesin]
        
        data$type = factor(data$type, levels = substitutions)
        
        
        if (nrow(data)==0)
        {
          warning(sprintf("No variants found for mutation type %s", mut))
          next
        }
        
        # These variables will be available at run-time, but not at compile-time.
        # To avoid compiling trouble, we initialize them to NULL.
        location = NULL
        
        # make rainfall plot
        plot = ggplot(data, aes(x=location, y=distance)) +
          geom_point(aes(colour=factor(type)), cex=cex) + 
          geom_vline(xintercept = as.vector(chr_cum), linetype="dotted") +
          annotate("text", x = m, y = ylim, label = labels, cex=cex_text) +
          xlab("Genomic Location") +
          ylab("Genomic Distance") +
          scale_y_log10() +
          scale_colour_manual(values=colors) +
          scale_x_continuous(expand = c(0,0), limits=c(0, max(chr_cum))) +
          ggtitle(title) +
          theme_bw() +
          theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) + 
          guides(colour = guide_legend(nrow = 2),
                 fill = guide_legend(nrow = 2))
        
        plots[[mut]] = plot
      }
      
      return(plot_grid(plotlist=plots, align="v", ncol = 1))
    } else if (method == "combine")
    {
      muts = unname(unlist(muts))
      vcf = vcf[which(1:length(vcf) %in% muts),]
  
      colors = unname(unlist(colors_list[type]))
  
      # get chromosome lengths of reference genome
      chr_length = seqlengths(vcf)
  
      # subset
      chr_length = chr_length[names(chr_length) %in% chromosomes]
  
      # cumulative sum of chromosome lengths
      chr_cum = c(0, cumsum(as.numeric(chr_length)))
  
      # Plot chromosome labels without "chr"
      names(chr_cum) = names(chr_length)
      labels = gsub("chr", "", names(chr_length))
  
      # position of chromosome labels
      m=c()
      for(i in 2:length(chr_cum))
          m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)
  
      # mutation characteristics
      type = loc = dist = chrom = c()
  
      # for each chromosome get types, location, distance and chromosome information
      # of all mutations
      for(i in 1:length(chromosomes))
      {
          type_list = list()
          chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
          
          n = length(chr_subset)
          if(n<=1){next}
          if ("snv" %in% type)
            type_list[["snv"]] = mut_type(chr_subset, type = type)
          if ("dbs" %in% type)
          {
            if (!isEmpty(type_list[["dbs"]]))
              for (j in 1:length(type_list[["dbs"]])) type_list[["dbs"]][j] = paste0(substr(type_list[["dbs"]][j],1,2),">NN")
          } else if ("indel" %in% type)
          {
            type_list[["indel"]] = mut_context(chr_subset, ref_genome, type = mut)
            type_list[["indel"]] = do.call(rbind, strsplit(type_list[["indel"]], "\\."))[,c(1,2,4)]
            type_list[["indel"]] = list("indel"=paste(type_list[["indel"]][,1], type_list[["indel"]][,2], type_list[["indel"]][,3], sep="."))
          }
          type = c(type, unname(unlist(type_list))[-1])
          loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
          dist = c(dist, diff(start(chr_subset)))
          chrom = c(chrom, rep(chromosomes[i],n-1))
      }
  
      data = data.frame(type = type,
                          location = loc,
                          distance = dist,
                          chromosome = chrom)
  
      # Removes colors based on missing mutation types.  This prevents colors from
      # shifting when comparing samples with low mutation counts.
      substitutions = unname(unlist(substitutions_list[mut_type]))
      typesin = substitutions %in% levels(data$type)
      colors = colors[typesin]
      
      data$type = factor(data$type, levels = substitutions)
  
      # These variables will be available at run-time, but not at compile-time.
      # To avoid compiling trouble, we initialize them to NULL.
      location = NULL
  
      # make rainfall plot
      plot = ggplot(data, aes(x=location, y=distance)) +
          geom_point(aes(colour=factor(type)), cex=cex) + 
          geom_vline(xintercept = as.vector(chr_cum), linetype="dotted") +
          annotate("text", x = m, y = ylim, label = labels, cex=cex_text) +
          xlab("Genomic Location") +
          ylab("Genomic Distance") +
          scale_y_log10() +
          scale_colour_manual(values=colors) +
          scale_x_continuous(expand = c(0,0), limits=c(0, max(chr_cum))) +
          ggtitle(title) +
          theme_bw() +
          theme(
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.key = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) + 
          guides(colour = guide_legend(nrow = 1))
  
      return(plot)
    }
}
