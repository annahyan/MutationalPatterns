#' Read VCF files into a GRangesList
#'
#' This function reads Variant Call Format (VCF) files into a GRanges object
#' and combines them in a GRangesList.  In addition to loading the files, this
#' function applies the same seqlevel style to the GRanges objects as the
#' reference genome passed in the 'genome' parameter.
#'
#' @param vcf_files Character vector of VCF file names
#' @param sample_names Character vector of sample names
#' @param genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome of your VCFs
#' @param group (Optional) Selector for a seqlevel group.  All seqlevels outside
#'              of this group will be removed.  Possible values:
#'              * 'all' for all chromosomes;
#'              * 'auto' for autosomal chromosomes;
#'              * 'sex' for sex chromosomes;
#'              * 'auto+sex' for autosomal + sex chromosomes (default);
#'              * 'circular' for circular chromosomes;
#'              * 'none' for no filtering, which results in keeping all
#'                seqlevels from the VCF file.
#' @param check_alleles (Optional) logical. If TRUE (default) positions with insertions,
#'              deletions and/or multiple alternative alleles are excluded
#'              from the vcf object, since these positions cannot be analysed
#'              with this package.  This setting can be set to FALSE to speed
#'              up processing time only if the input vcf does not contain any
#'              of such positions, as these will cause obscure errors.
#' @param n_cores (Optional) numeric. Number of cores used for parallel processing. If no
#'              value is given, then the number of available cores is autodetected.
#' 
#' @return A GRangesList containing the GRanges obtained from 'vcf_files'
#'
#' @importFrom BiocGenerics match
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb "seqlevelsStyle<-"
#' @importFrom GenomeInfoDb "organism"
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom GenomeInfoDb extractSeqlevelsByGroup
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @importFrom plyr llply
#'
#' @examples
#' # The example data set consists of three colon samples, three intestine
#' # samples and three liver samples.  So, to map each file to its appropriate
#' # sample name, we create a vector containing the sample names:
#' sample_names <- c ( "colon1", "colon2", "colon3",
#'                     "intestine1", "intestine2", "intestine3",
#'                     "liver1", "liver2", "liver3" )
#'
#' # We assemble a list of files we want to load.  These files match the
#' # sample names defined above.
#' vcf_files <- list.files(system.file("extdata", 
#'                                     package="MutationalPatterns"),
#'                                     pattern = ".vcf", full.names = TRUE)
#'
#' # Get a reference genome BSgenome object.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library("BSgenome")
#' library(ref_genome, character.only = TRUE)
#'
#' # This function loads the files as GRanges objects
#' vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
#'
#' @export

read_vcfs_as_granges <- function(vcf_files, sample_names, genome,
                                    group = "auto+sex", check_alleles = TRUE, n_cores)
{
    # Check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Please provide the same number of sample names as VCF files")

    ref_genome <- base::get(genome)
    ref_organism <- GenomeInfoDb::organism(ref_genome)
    ref_style <- seqlevelsStyle(ref_genome)

    # Name the VCF's genome as the name of the genome build instead of
    # the BSgenome package name.
    genome_name <- genome(ref_genome)[[1]]

    # Check the class of the reference genome
    if (!(class(ref_genome) == "BSgenome"))
        stop("Please provide the name of a BSgenome object.")

    # If number of cores is not provided, detect the number of available cores.  
    # Windows does not support forking, only threading, so unfortunately, we 
    # have to set it to 1.
    # On confined OS environments, this value can be NA, and in such
    # situations we need to fallback to 1 core.

    if(missing(n_cores))
    {
      n_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(n_cores)))
          n_cores <- detectCores()
        else
          n_cores = 1
    }

    # We handle errors and warnings separately for mclapply, because the error
    # reporting of mclapply is done through its return value(s).
    original_warn_state = getOption("warn")
    options(warn=-1)

    # Store warning messages in a vector.
    warnings <- NULL

    # Show the warning once for all VCF files that are loaded with this
    # call to read_vcfs_as_granges.
    if (!check_alleles)
    {
        warning(paste("check_alleles is set to FALSE.  Make sure your",
                      "input VCF does not contain any positions with",
                      "multiple alternative",
                      "alleles, as these positions cannot be analysed",
                      "with MutationalPatterns and cause obscure",
                      "errors."), immediate. = TRUE)
    }

    vcf_list <- mclapply (seq_along(vcf_files), function (index)
    {
        file <- vcf_files[index]

        # Use VariantAnnotation's readVcf, but only store the
        # GRanges information in memory.  This speeds up the
        # loading significantly.
        vcf <- rowRanges(readVcf (file, genome_name))
        
        if (length(vcf) == 0)
          stop(sprintf("Vcf file %s is empty", file))
        
        # Convert to a single naming standard.
        seqlevelsStyle(vcf) <- ref_style[1]

        groups <- c()
        if (group != "none")
        {
            if (group == "auto+sex")
            {
                groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                                    style = ref_style,
                                                    group = "auto"),
                            extractSeqlevelsByGroup(species = ref_organism,
                                                    style = ref_style,
                                                    group = "sex"))

                # In some cases, the seqlevelsStyle returns multiple styles.
                # In this case, we need to do a little more work to extract
                # a vector of seqlevels from it.
                groups_names <- names(groups)
                if (! is.null(groups_names))
                {
                    # The seqlevels in the groups are now duplicated.
                    # The following code deduplicates the list items, so that
                    # creating a data frame will work as expected.
                    unique_names <- unique(groups_names)
                    groups <- llply(unique_names, function(x) groups[groups_names == x])
                    groups <- llply(groups, unlist, recursive = FALSE)

                    # In case there are multiple styles applied, we only use the first.
                    groups <- unique(as.vector(groups[[1]]))
                }
            }
            else
            {
                groups <- extractSeqlevelsByGroup ( species = ref_organism,
                                                   style = ref_style,
                                                   group = group )
                groups <- unique(as.vector(t(groups)))
            }

            # The provided VCF files may not contain all chromosomes that are
            # available in the reference genome.  Therefore, we only take the
            # chromosomes that are actually available in the VCF file,
            # belonging to the filter group.
            groups <- intersect(groups, seqlevels(vcf))

            # We use 'pruning.mode = "tidy"' to minimize the deleterious effect
            # on variants, yet, remove all variants that aren't in the filter
            # group.  By default, keepSeqlevels would produce an error.
            vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
        }

        if (check_alleles)
        {
            # Find and exclude positions with indels or multiple
            # alternative alleles.
            rem <- which(lengths(vcf$ALT) > 1)

            if (length(rem) > 0)
            {
                vcf = vcf[-rem]
                warnings$check_allele <- 
                  rbind(warnings$check_allele,
                        c(sample_names[[index]], length(rem)))
            }
        }
        
        # Search for DBS which are given as two sequential locations
        dbs = intersect(which(diff(start(vcf)) == 1),
                        which(nchar(as.character(vcf$REF)) == 1 &
                                nchar(as.character(unlist(vcf$ALT))) ==1 ))
        if (length(dbs) > 0)
          if (dbs[length(dbs)] == length(vcf))
            dbs <- dbs[-length(dbs)]
          mnv = NULL
          if (1 %in% diff(dbs))
          {
            # Remove multi nucleotide variants
            rem = which(diff(dbs) == 1)
            mnv = unique(sort(c(dbs[rem],dbs[rem]+1,dbs[rem]+2)))
            
            rem <- unique(c(rem,rem+1))
            dbs <- dbs[-rem]
            warnings$dbs <- rbind(warnings$dbs,
                                  c(sample_names[[index]], length(mnv)))
          }
        
        # If there are DBS, then change end position of variant,
        # add second ref and second alt base and delete next variant from vcf
        if (length(dbs) > 0)
        {
          end(vcf)[dbs] = start(vcf)[dbs]+1
          vcf$REF[dbs] = DNAStringSet(paste0(as.character(vcf$REF[dbs]), as.character(vcf$REF[dbs+1])))
          vcf$ALT[dbs] = DNAStringSetList(lapply(dbs, function(i) {
            DNAStringSet(paste0(as.character(unlist(vcf$ALT[i])), as.character(unlist(vcf$ALT[i+1]))))
            }))
          vcf = vcf[-c(mnv,(dbs+1)),]
        }
        
        indel = which((nchar(as.character(vcf$REF)) == 1 &
                        nchar(as.character(unlist(vcf$ALT))) > 1) |
                        (nchar(as.character(vcf$REF)) > 1 &
                           nchar(as.character(unlist(vcf$ALT))) == 1))
          
        if (any(c(grepl("[^ACGT]",as.character(vcf$REF[indel])),
                  grepl("[^ACGT]",as.character(unlist(vcf$ALT[indel])))))){
          vcf = vcf[-indel,]
          warnings$indel <- rbind(warnings$indel,
                                  c(sample_names[[index]]))
        }

        # Pack GRanges object and the warnings to be able to display warnings
        # at a later time.
        return(list(vcf, warnings))
    }, mc.cores = n_cores)

    # Reset the option.
    options(warn=original_warn_state)

    # mclapply wraps the call into a try(..., silent=TRUE)
    # When an error occurs, the error is returned, and accessible in the
    # return value(s).  The for-loop below checks for erroneous returns
    # and shows the error message of the first occurring error.
    #
    # The return values of the mclapply output are packed as
    # list(<GRanges>, <warnings>).  The function below unpacks the GRanges
    # and displays the warnings.
    
    # Handle errors
    summ <- lapply(vcf_list, function(item) {
        if (class(item) == "try-error") stop (item)
        ref = as.character(item[[1]]$REF)
        alt = as.character(unlist(item[[1]]$ALT))
        
        nsnvs = length(which(nchar(ref) == 1 & nchar(alt) == 1))
        ndbs = length(which(nchar(ref) == 2 & nchar(alt) == 2))
        nindel = length(which((nchar(ref) == 1 & nchar(alt) != 1) | 
                                (nchar(ref) != 1 & nchar(alt) == 1)))
        
        return(c(nsnvs, ndbs, nindel))
    })
    
    
    # Handle warnings
    warnings <- do.call(rbind, vcf_list)[,2]
    warnings <- sapply(warnings, function(item) {
      if (is.null(item))
        item = list("check_alleles"=NULL,"dbs"=NULL,"indel"=NULL)
      item[c('check_alleles','dbs','indel')]
      })
    warns = NULL
    for (i in which(!(is.na(rownames(warnings))))){
        warns[[rownames(warnings)[i]]] = do.call(rbind, warnings[i,])
        colnames(warns[[rownames(warnings)[i]]]) <- c("Sample", "Position(s)")
    }
  
    lapply(names(warns), function(item) {
        if (item == "check_alleles")
        {
            warning("Position(s) with multiple ",
                    "alternative alleles are excluded\n",
                    paste0(capture.output(warns$check_allele), collapse = "\n"),
                    immediate. = TRUE)
        } else if (item == "dbs")
        {
          warning("Position(s) excluded that form ", 
                  "multiple nucleotide variants ", 
                  "of length more than 2.\n",
                  paste0(capture.output(warns$dbs), collapse = "\n"),
                  immediate. = TRUE)
        } else if (item == "indel")
        {
          warning("Indels not according to VCF format ",
                  "4.2 or higher, all indels are excluded from:\n ",
                  paste(warns$indel, collapse=", "),
                  immediate. = TRUE)
        }
    })
    
    vcf_list <- GRangesList(do.call(rbind, vcf_list)[,1])
    
    # Set the provided names for the samples.
    names(vcf_list) <- sample_names
    
    # Print a summary of mutations found
    summ <- do.call(rbind, summ)
    colnames(summ) <- c("Number of SNV", "Number of DBS", "Number of indel")
    rownames(summ) <- names(vcf_list)
    print(summ)
    
    return(vcf_list)
}
