#' Filter OTU table by occurrence and abundance
#'
#' @param ps
#' @param  occurrence.threshold
#' @param  abundance.threshold
#'
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export


filter_ps  = function(ps = ps,  occurrence.threshold=0 ,  abundance.threshold= 0 ){
ps1<- filter_taxa(ps,function(x)sum(x>0)>occurrence.threshold,TRUE)
ps2 <- prune_taxa(taxa_sums(ps1)>abundance.threshold,ps1)
#filteredps<-rarefy_even_depth(ps2,  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
filteredps=ps2
return(filteredps)
}
