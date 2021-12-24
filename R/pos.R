#' Create function that averages only positive values in a vector
#' @param  vector
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}
