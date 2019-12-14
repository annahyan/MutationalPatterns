#' Count occurrences per base substitution type and strand
#' 
#' For each base substitution type and strand the total number
#' of mutations and the relative contribution within a group is returned.
#'
#' @param mut_mat_s 192 feature mutation count matrix, result from
#' 'mut_matrix_stranded()'
#' @param by Character vector with grouping info, optional
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

strand_occurrences = function(mut_mat_s, by)
{
    df = t(mut_mat_s)

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
    strand = unlist(lapply( strsplit(rownames(mut_mat_s), "-") , function(x) x[2]))
    # get substitutions from rownames
    substitutions = substring(rownames(mut_mat_s), 3, 5)
    
    # sum per substition per strand
    x2 = melt(aggregate(t(x), by = list(substitutions, strand), FUN=sum))
    x2_r = melt(aggregate(t(x_r), by = list(substitutions, strand), FUN=sum))
    colnames(x2) = c("type", "strand", "group", "no_mutations")
    colnames(x2_r) = c("type", "strand", "group", "relative_contribution")

    # combine relative and absolute
    y = merge(x2, x2_r)

    # reorder group, type, strand
    y = y[,c(3,1,2,4,5)]
    y = y[order(y$group, y$type),]

    return(y)
}
