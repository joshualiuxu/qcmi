#' Find the number of zeroes in a vector
#' @param vec
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}
