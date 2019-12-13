#' Retrieve context of indel
#'
#' A function to extract the indel context from a data.frame, 
#' specified by the context according to the COSMIC database
#' @param indel A data frame with information about a indel. Columns 
#' needed in data.frame are: \cr
#' \itemize{
#' \item{indel_len:} {  Length of an indel sequence}
#' \item{indel_seq:} {  Sequence of an indel}
#' \item{indel_type:} {  'del'(etion) or 'ins'(erstion)}
#' \item{repeats:} {  Number of repeats of an indel}
#' \item{bimh:} {  Number of bases in microhomology}
#' }
#' @return Character vector with the context of the indel
#' @import Biostrings
#'
#' @examples
#' ## Get a data.frame with the columns needed. Is also part of the
#' ## function extract_indels
#' indel = data.frame("chrom" = "chr1",
#'                    "pos" = 3249117,
#'                    "ref" = "AT",
#'                    "alt" = "A",
#'                    "indel_len" = 1,
#'                    "indel_type" = "del",
#'                    "indel_seq" = "T",
#'                    "repeats" = 4,
#'                    "bimh" = 1)
#' 
#' context = cosmic_indel_context(indel)
#'
#' @seealso
#' \code{\link{extract_indels}}
#'
#' @export

cosmic_indel_context = function(indel) 
{
    len = indel$indel_len
    sq = as.character(indel$indel_seq)
    type = indel$indel_type
    n = indel$repeats
    bimh = indel$bimh
    
    # Get C or T base of indel
    if (sq == "A" | sq == "G") 
        sq = as.character(complement(DNAString(sq)))
    
    # If type is a deletion, look at number of repeats
    # and look at number of bases in microhomology
    if (type == "del")
    {
        if (n >= 2)
        {
            if (len == 1) 
            { 
              len = paste0("1bp.homopol.",sq,".len") 
            } else if (len >= 5) 
            { 
              len = "rep.len.5+.rep" 
            } else 
            { 
              len = paste0("rep.len.",len,".rep") 
            }
            
            if (n >= 6) { n = "6+" }
            else { n = as.character(n) }
          
            return(sprintf("del.%s.%s", len, n))
        } else {
            if (len == 1)
            {
              return(sprintf("del.1bp.homopol.%s.len.1", sq))
            }
            if (len >= 5) { len = "5+"}
          
            if (bimh >= 1)
            {
                if (bimh >= 5){ bimh = "5+" }
                else { bimh = as.character(bimh) }
                
                return(sprintf("del.mh.len.%s.bimh.%s", len, bimh))
            } else {
                return(sprintf("del.rep.len.%s.rep.1", len))
            }
        }
    # If type is insertion, look at number of repeats. 
    # No microhomology context for insertions
    } else if (type == "ins") 
    { 
        n = n-1
    
        if (n >= 1)
        {
            if (len == 1) 
            { 
                len = paste0("1bp.homopol.",sq,".len") 
            } else if (len >= 5) 
            { 
                len = "rep.len.5+.rep" 
            } else 
            { 
                len = paste0("rep.len.",len,".rep") 
            }
            
            if (n >= 5) { n = "5+" }
            else { n = as.character(n) }
        
            return(sprintf("ins.%s.%s", len, n))
        } else {
            if (len == 1)
            {
                return(sprintf("ins.1bp.homopol.%s.len.0", sq))
            } else if (len >= 5) 
            {
                return("ins.rep.len.5+.rep.0")
            } else {
                return(sprintf("ins.rep.len.%s.rep.0", len))
            }
        }
    }
}
