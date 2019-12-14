#' Test for enrichment or depletion of mutations in genomic regions
#'
#' This function aggregates mutations per group (optional) and performs an
#' enrichment depletion test.
#'
#' @param x data.frame result from genomic_distribution() 
#' @param by Optional grouping variable, e.g. tissue type
#' @return data.frame with the observed and expected number of mutations per
#' genomic region per group (by) or sample
#'
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#' @importFrom stats aggregate
#'
#' @examples
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' distr <- readRDS(system.file("states/distr_data.rds",
#'                     package="MutationalPatterns"))
#' 
#' tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#'
#' ## Perform the enrichment/depletion test by tissue type.
#' distr_test <- enrichment_depletion_test(distr, by = tissue)
#'
#' ## Or without specifying the 'by' parameter.
#' distr_test2 <- enrichment_depletion_test(distr)
#'
#' @seealso
#' \code{\link{genomic_distribution}},
#' \code{\link{plot_enrichment_depletion}}
#'
#' @export

enrichment_depletion_test = function(x, by = c())
{
    # Handle the 'by' parameter when necessary by aggregating x
    if (length(by) > 0)
    {
        x$by = by
        # Sum the columns while aggregating rows based on unique values
        # in 'by' and 'region'.
        res2 = stats::aggregate(cbind(n_muts,
                                        surveyed_length,
                                        surveyed_region_length,
                                        observed) ~ by + region,
                                data = x, sum)
    }
    else
    {
        res2 = x
        # In this case, the 'by' variable is 'sample' variable.
        res2$by = res2$sample
        # Select output columns
        res2 = res2[,c(9,1,3,4,6,8)]
    }

    # Calculate probability and expected number of mutations
    res2$prob = res2$n_muts / res2$surveyed_length
    res2$expected = res2$prob * res2$surveyed_region_length

    # Perform enrichment/depletion test for each row
    res3 = data.frame()
    for(i in 1:nrow(res2))
    {
        x = res2[i,]
        res3 = rbind(res3, binomial_test(x$prob,
                                            x$surveyed_region_length,
                                            x$observed))
    }

    # Combine results into one data frame
    df = cbind(res2, res3)
    return(df)
}
