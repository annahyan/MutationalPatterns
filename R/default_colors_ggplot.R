#' Default colors from ggplot2 package
#' 
#' Get the default colors in HEX format which are used by the ggplot2 
#' package, such that the colors can be used for the same signatures, 
#' even if the order or number of signatures differs. \cr\cr
#' Colors are chosen by equal distribution of the range [0,360]
#' over the number of colors wanted. 0 yields red, 120 yields green,
#' 240 yields blue etc.
#' 
#' @param n Numeric stating the number of default colors wanted
#' @return Character vector with colornames in HEX format
#' 
#' @importFrom grDevices hcl
#' 
#' @examples 
#' ## Get the default colors of ggplot2 for 10 classes
#' n = 10
#' 
#' default_colors = default_colors_ggplot(n)
#' 
#' @export

default_colors_ggplot <- function(n) {
  
  # Pick numbers within [0,360] with equal distribution
  hues = seq(15, 375, length = n + 1)
  
  # Translate angles in HEX format
  return(hcl(h = hues, l = 65, c = 100)[1:n])

}

