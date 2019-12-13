#' Count occurrences per base substitution type and strand
#' 
#' For each base substitution type and strand the total number
#' of mutations and the relative contribution within a group is returned.
#'
#' @param mut_mat_s strand feature mutation count matrix, result from
#' 'mut_matrix_stranded()'
#' @param by (Optional) Character vector with grouping info
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param method (Optional) Character stating if all mutation types should be combined
#' (='combine') in relative contribution or taking apart (= 'split'). 
#' Default is 'split'
#'
#' @return A data.frame with the total number of mutations and relative
#' contribution within group per base substitution type and strand 
#'
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
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
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{plot_strand}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_occurrences = function(mut_mat_s, by, type, method = "split")
{
    # Check the mutation type argument
    type = check_mutation_type(type)
    
    # If mutation matrix object is a matrix, then find the mutation type
    if (class(mut_mat_s) == "matrix")
    {
      if ((all(rownames(mut_mat_s) %in% TRIPLETS_192_trans) | 
          all(rownames(mut_mat_s) %in% TRIPLETS_192_rep)) &
          length(rownames(mut_mat_s) > 0))
      {
        mut_mat_s = list("snv" = mut_mat_s)
        type = "snv"
      } else if ((all(rownames(mut_mat_s) %in% DBS_trans) |
                all(rownames(mut_mat_s) %in% DBS_rep)) &
               length(rownames(mut_mat_s)) > 0)
      {
        mut_mat_s = list("dbs" = mut_mat_s)
        type = "dbs"
      } else if (all(unique(type_context$context) %in% INDEL_CONTEXT))
      {
        type_context = list("indel"=list("types"=type_context$types,
                                         "context"=type_context$context))
        type = "indel"
      }
    }
    
    # Get the asked types
    type = intersect(type, names(mut_mat_s))
  
    if (method == "split")
    {
      z = list()
      
      for (m in type)
      {
        df = t(mut_mat_s[[m]])
        
        # check if grouping parameter by was provided, if not group by all
        if(missing(by)){by = rep("all", nrow(df))}
        
        # sum by group
        x = stats::aggregate(df, by=list(by), FUN=sum) 
        
        # add group as rownames
        rownames(x) = x[,1]
        x = x[,-1]
        
        # calculate relative contribution within group
        x_r = x / rowSums(x)
        
        # get strand from rownames
        strand = unlist(lapply( strsplit(rownames(mut_mat_s[[m]]), "-") , function(x) x[2]))
        # get substitutions from rownames
        if (m == "snv"){ substitutions = substring(rownames(mut_mat_s[[m]]), 3, 5) }
        else if (m == "dbs"){ substitutions = substring(rownames(mut_mat_s[[m]]), 1, 2) }
        else if (m == "indel") 
        { 
          subs = do.call(rbind, strsplit(rownames(mut_mat_s[[m]]), "\\."))[,1:4]
          substitutions = paste(subs[,1],subs[,2],subs[,3],subs[,4],sep=".")
        }
        
        # sum per substition per strand
        x2 = melt(aggregate(t(x), by = list(substitutions, strand), FUN=sum))
        x2_r = melt(aggregate(t(x_r), by = list(substitutions, strand), FUN=sum))
        colnames(x2) = c("type", "strand", "group", "no_mutations")
        colnames(x2_r) = c("type", "strand", "group", "relative_contribution")
        
        # combine relative and absolute
        y = merge(x2, x2_r)
        
        # reorder group, type, strand
        y$mutation = m
        y = y[,c(3,6,1,2,4,5)]
        y = y[order(y$group, y$type),]
        
        z[[m]] = y
      }
      
      # Return a vector when there is only 1 mutation type
      if (length(names(z)) == 1)
        z = z[[1]]
      
      return(z)
    } else if (method == "combine"){
      mutation_types = names(mut_mat_s)
      df = t(do.call(rbind, mut_mat_s))
      
      # check if grouping parameter by was provided, if not group by all
      if(missing(by)){by = rep("all", nrow(df))}
      
      # sum by group
      x = stats::aggregate(df, by=list(by), FUN=sum) 
      
      # add group as rownames
      rownames(x) = x[,1]
      x = x[,-1]
      
      # calculate relative contribution within group
      x_r = x / rowSums(x)
      
      # get strand from rownames
      strand = unlist(lapply( strsplit(colnames(df), "-") , function(x) x[2]))
      # get substitutions from rownames
      substitutions = c()
      mut_types = c()
      for (m in mutation_types)
      {
        if (m == "snv"){ 
          substitutions = c(substitutions, substring(rownames(mut_mat_s[[m]]), 3, 5))
          mut_types = c(mut_types, rep("snv", nrow(mut_mat_s[[m]])))
        } else if (m == "dbs")
        { substitutions = c(substitutions, substring(rownames(mut_mat_s[[m]]), 1, 2))
          mut_types = c(mut_types, rep("dbs", nrow(mut_mat_s[[m]])))
        }
      }
      
      # sum per substition per strand
      x2 = melt(aggregate(t(x), by = list(mut_types, substitutions, strand), FUN=sum))
      x2_r = melt(aggregate(t(x_r), by = list(mut_types, substitutions, strand), FUN=sum))
      colnames(x2) = c("mutation", "type", "strand", "group", "no_mutations")
      colnames(x2_r) = c("mutation", "type", "strand", "group", "relative_contribution")
      
      # combine relative and absolute
      y = merge(x2, x2_r)
      
      # reorder group, type, strand
      y = y[,c(4,1,2,3,5,6)]
      y = y[order(y$mutation, decreasing = TRUE),]
      
      return(y)
    } else {stop("Unknown value for 'method' is given. Choose between 'split' (default) or 'combine'")}
}
