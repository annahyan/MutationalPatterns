#' Significance test for strand asymmetry
#'
#' This function performs a Poisson test for the ratio between mutations on 
#' each strand
#' 
#' @param strand_occurrences Dataframe with mutation count per strand, result
#' from strand_occurrences()
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param method (Optional) Character stating how to use the data. 
#' \itemize{
#'   \item{"split":} { Each mutation type has a seperate strand bias test}
#'   \item{"combine":} { Combined strand bias test for all mutation types}
#' }   
#' Default is "split"
#' @return Dataframe with poisson test P value for the ratio between the
#' two strands per group
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom plyr .
#' @importFrom plyr ddply
#' @importFrom plyr summarise
#' @importFrom stats "poisson.test"
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
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_bias_test = function(strand_occurrences, type, method = "split")
{
    # Check mutation type argument
    type = check_mutation_type(type)
    
    if(class(strand_occurrences) == "data.frame")
    {
      if(length(unique(strand_occurrences$mutation)) > 1)
        warning(paste("No named list found for 'strand_occurrences'.",
                      "Method is set to 'combine'"))
      method = "combine"        
    }
    
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    group = NULL
    strand = NULL
    variable = NULL

    # statistical test for strand ratio
    # poisson test
    
    if (method == "split")
    {
      # Get the asked mutation types
      type = intersect(type, names(strand_occurrences))
      
      df_result = list()
      
      # For each type, perform the strand bias test
      for (m in type)
      {
        df_strand = reshape2::dcast(melt(strand_occurrences[[m]]),
                                    group + mutation + type ~ strand,
                                    sum,
                                    subset = plyr::.(variable == "no_mutations"))
        
        df_strand$total = df_strand[,4] + df_strand[,5]
        df_strand$ratio = df_strand[,4] / df_strand[,5]
        df_strand$p_poisson = apply(df_strand, 1, function(x) poisson.test(c(as.numeric(x[4]), as.numeric(x[5])), r=1)$p.value)
        df_strand$significant[df_strand$p_poisson < 0.05] = "*"
        df_strand$significant[df_strand$p_poisson >= 0.05] = " "
        df_result[[m]] = df_strand
      }
    } else if (method == "combine")
    {
      # Perform the strand bias test
      df_strand = reshape2::dcast(melt(strand_occurrences),
                                  group + mutation + type ~ strand,
                                  sum,
                                  subset = plyr::.(variable == "no_mutations"))
      
      df_strand$total = df_strand[,4] + df_strand[,5]
      df_strand$ratio = df_strand[,4] / df_strand[,5]
      df_strand$p_poisson = apply(df_strand, 1, function(x) poisson.test(c(as.numeric(x[4]), as.numeric(x[5])), r=1)$p.value)
      df_strand$significant[df_strand$p_poisson < 0.05] = "*"
      df_strand$significant[df_strand$p_poisson >= 0.05] = " "
      df_result = df_strand
    }

    return(df_result)
}
