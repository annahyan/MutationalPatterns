#' Set global context, class and color values for indels
#' 
#' A function to set the global values of the context variable, 
#' the class variable and the color variable of indels. For indels
#' there is no single intuitive and naturally constrained set of 
#' mutation types, so a custom set of mutation types can be used
#'
#' @param indel (Optional) List of vectors for context, class and color of indels.
#' Color vector must have the same length as the class vector, 
#' because in plotting profiles, each class is represented by one color.\cr\cr
#' It is also possible to give a character for predefined variables:
#' \itemize{
#'   \item{"predefined"} {Represents a indel context of 3 classes per 
#'   deletion and indel: "repetitive region", "microhomology" and 
#'   "none" of these two. Indels have lengths 1 to 5+}
#'   \item{"cosmic"} {Represents the indel context according to the
#'   COSMIC database}
#' }
#' Default is "cosmic"
#' @return This function changes internal variables used by the package
#' to analyse indels.
#' 
#'
#' @examples 
#' ## For custom input of predefined context
#' indel = list("matrix" = mut_matrix,
#'              "context" =c(
#' paste0('del.rep.len.', 1:5),
#' paste0('ins.rep.len.', 1:5),
#' paste0('del.mh.bimh.', 1:5),
#' paste0('ins.mh.bimh.', 1:5),
#' paste0('del.none.len.', 1:5),
#' paste0('ins.none.len.', 1:5)
#' ),
#'              "class"=c(
#' rep('del.rep', 5), rep('ins.rep', 5),
#' rep('del.mh', 5), rep('ins.mh', 5),
#' rep('del.none', 5), rep('ins.none', 5)
#' ),
#'              "colors" = c(
#' "#F7BF80", "#ED8212", 
#' "#B5D988", "#31A12C", 
#' "#E44A39", "#B81C20"
#' ))
#' 
#' indel_mutation_type(indel)
#'
#' @export

indel_mutation_type <- function(indel)
{
  # Default indel context is "cosmic"
  if (missing(indel)) { indel = "cosmic" }
  
  # If a character is given, check for "predefined" or "cosmic"
  # and set global variables accordingly

  e <- loadNamespace("MutationalPatterns")
  unlockBinding("INDEL",e)
  unlockBinding("INDEL_CONTEXT",e)
  unlockBinding("INDEL_CLASS",e)
  unlockBinding("INDEL_CLASS_HEADER",e)
  unlockBinding("COLORS_INDEL",e)
  
  if (class(indel) == "character")
  {
    if (indel == "cosmic"){
      if (e$INDEL != "cosmic")
      {
        e$INDEL = "cosmic"
        e$INDEL_CONTEXT = e$INDEL_CONTEXT_COSMIC
        e$INDEL_CLASS = e$INDEL_CLASS_COSMIC
        e$INDEL_CLASS_HEADER =  e$INDEL_CLASS_HEADER_COSMIC
        e$COLORS_INDEL = e$COLORS_INDEL_COSMIC
      }
    } else if (indel == "predefined"){
      e$INDEL = "predefined"
      e$INDEL_CONTEXT = e$INDEL_CONTEXT_PREDEF
      e$INDEL_CLASS = e$INDEL_CLASS_PREDEF
      e$INDEL_CLASS_HEADER =  e$INDEL_CLASS_HEADER_PREDEF
      e$COLORS_INDEL = e$COLORS_INDEL_PREDEF
    } else { stop(print("Unknown option for indel mutation types",
                        "Provide either 'predefined' or 'cosmic' character or",
                        "list with indel classes and contexts")) }
  } else if (class(indel) == "list")
  {
    # If a list is given, then obligatory elements of list are 
    # the mutation count matrix, the classes and contexts of the indels
    needed = c("matrix", "class","context")
    if (all(needed %in% names(indel)))
    {
      if ("header" %in% names(indel)) { e$INDEL_CLASS_HEADER = indel$header }
      else { e$INDEL_CLASS_HEADER = NULL }
      
      # Colors can be given in the list, else default colors are chosen
      if ("colors" %in% names(indel)) { e$COLORS_INDEL = indel$colors }
      else 
      {
        e$COLORS_INDEL = default_colors_ggplot(length(unique(indel$class)))
      }
      e$INDEL_MATRIX = indel$matrix
      e$INDEL_CONTEXT = indel$context
      e$INDEL_CLASS = indel$class
      e$INDEL = "custom"
    } else 
    {
      # Stop when one or more obligatory elements are missing
      stop(sprintf("Custom indel list misses information about: %s", 
                   needed[which(!(needed %in% names(indel)))] ))
    }
  }
  
  lockBinding("INDEL",e)
  lockBinding("INDEL_CONTEXT",e)
  lockBinding("INDEL_CLASS",e)
  lockBinding("INDEL_CLASS_HEADER",e)
  lockBinding("COLORS_INDEL",e)
}
