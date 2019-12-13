#' Extract indel sequence, type and length
#'
#' @param bed A dataframe containing the columns: chrom, pos, ref, alt
#' @param ref_genome (Optional) A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param get_other_indel_allele (Optional) Only applies when mode=='indel' For indels, some vcfs only report
#' the sequence of one allele (REF for deletions and ALT for insertions). If TRUE, the unreported
#' allele will be retrieved from the genome: a 5' base relative to the indel sequence. This base
#' will also be added to the indel sequence and the POS will be adjusted accordingly (POS=POS-1).
#' @param verbose (Optional) Print progress messages?
#'
#' @importFrom VariantAnnotation alt
#' @importFrom VariantAnnotation ref
#'
#' @return A dataframe in the same structure as a bed file
#' @export
get_contexts_indel <- function(bed, ref_genome=DEFAULT_GENOME, 
                               get_other_indel_allele=FALSE, verbose=FALSE){
  
  bed_colnames <- c('chrom','pos','ref','alt')
  if(!(identical(colnames(bed)[1:4], bed_colnames))){
    warning("colnames(bed)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
    colnames(bed)[1:4] <- bed_colnames
  }
  
  bed$chrom = as.character(bed$chrom)
  bed$ref = as.character(bed$ref)
  bed$alt = as.character(bed$alt)
  
  if(verbose){ message('Converting chrom name style to style in ref_genome...') }
  seqlevelsStyle(bed$chrom) <- seqlevelsStyle(eval(parse(text=ref_genome)))
  
  if(verbose){ message('Determining indel type...') }
  ## Calc sequence lengths
  bed$ref_len <- nchar(bed$ref)
  bed$alt_len <- nchar(bed$alt)
  
  ## Remove snvs
  bed <- bed[!(bed$ref_len==1 & bed$alt_len==1),]
  
  ## Determine indel type
  bed$indel_type <- with(bed,{
    unlist(Map(function(ref_len, alt_len){
      if(ref_len >= 2 & alt_len >= 2){
        if(ref_len == alt_len){ 'mnv_neutral' }
        else if(ref_len > alt_len){ 'mnv_del' }
        else if(ref_len < alt_len){ 'mnv_ins' }
      } else if(ref_len > alt_len){
        'del'
      } else if(ref_len < alt_len){
        'ins'
      }
    },ref_len, alt_len, USE.NAMES=FALSE))
  })
  
  if(get_other_indel_allele==TRUE){
    if(verbose){ message('Retrieving other indel allele...') }
    bed_split <- lapply(
      list(del_type=c('del','mnv_del'),ins_type=c('ins','mnv_ins'),mnv_neutral='mnv_neutral'),
      function(i){ bed[bed$indel_type %in% i, ] }
    )
    
    if(nrow(bed_split$del_type)!=0){
      ## Deletions
      ## ref:   'AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
      ## alt:  'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
      ## nchar('AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC') = nchar(ref) = 46
      ## getSeq(x=eval(parse(text = ref_genome)), names='chr4',start=84726292-1,84726292+46-1)
      ##       'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
      ## 5' base relative to ref ->
      ##    alt column
      ##    ref sequence
      bed_split$del_type$alt <- with(bed_split$del_type, {
        getSeq(
          x=eval(parse(text=ref_genome)),
          names=chrom, start=pos-1,end=pos-1,
          as.character=TRUE
        )
      })
      bed_split$del_type$ref <- with(bed_split$del_type, { paste0(alt, ref) })
      bed_split$del_type$pos <- bed_split$del_type$pos-1
    }
    
    if(nrow(bed_split$ins_type)!=0){
      ## Insertions
      ## ref:  ''
      ## alt:  'AGAGAGAGAGACAGAA'
      ## nchar('AGAGAGAGAGACAGAA') = nchar(alt) = 16
      ## getSeq(x=eval(parse(text = ref_genome)), names='chr12',start=6902128-1,6902128+16-1)
      ##      'GAGAGAGAGAGACAGAA'
      ## 5' base relative to alt ->
      ##    ref column
      ##    alt sequence
      ## Substract 1 from pos
      bed_split$ins_type$ref <- with(bed_split$ins_type, {
        getSeq(
          x=eval(parse(text=ref_genome)),
          names=chrom, start=pos-1,end=pos-1,
          as.character=TRUE
        )
      })
      bed_split$ins_type$alt <- with(bed_split$ins_type, { paste0(ref,alt) })
      bed_split$ins_type$pos <- bed_split$ins_type$pos-1
    }
    
    ## Unsplit bed
    bed <- do.call(rbind, bed_split)
    rownames(bed) <- NULL
    
    ## Recalculate ref/alt length
    bed$ref_len <- nchar(bed$ref)
    bed$alt_len <- nchar(bed$alt)
  }
  
  if(verbose){ message('Determining indel length and sequence...') }
  ## Determine indel length
  bed$indel_len <- abs(bed$alt_len-bed$ref_len)
  
  ## Determine indel seq
  bed$indel_seq <- with(bed,{
    unlist(Map(function(ref,alt,indel_type,indel_len){
      indel_start_pos <- 2
      if(indel_type %in% c('del','mnv_del')){ ## dels
        substring(ref, indel_start_pos, indel_start_pos+indel_len-1)
      } else if(indel_type %in% c('ins','mnv_ins')){ ## ins
        substring(alt, indel_start_pos, indel_start_pos+indel_len-1)
      } else {
        NA
      }
    },ref,alt,indel_type,indel_len, USE.NAMES=FALSE))
  })
  
  ## Output
  if(verbose){ message('Returning indel characteristics...') }
  out <- bed[bed$indel_type %in% c('ins','del'),]
  out <- out[,c('chrom','pos','ref','alt','indel_len','indel_type','indel_seq')]
  return(out)
}
