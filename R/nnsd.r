#' Accessory function for Random Matrix Theory
#' @param sp
#' @return 
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn} 
#' @references
#' @export


nnsd = function(sp){
			nns <- NULL
			for(j in 2:length(sp)){
				nn <- abs(sp[j] - sp[j-1])
				nns <- c(nns, nn)
			}
			nns
		}