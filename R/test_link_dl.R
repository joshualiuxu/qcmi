#' Classifying ecological associations of dispersal limitation
#' @param   link_table_row
#' @param   OTUabd
#' @param   geodist
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

test_link_dl = function(link_table_row, OTUabd, geodist){
  OTU1 = as.character(as.matrix(link_table_row[1]))
  OTU2	 = as.character(as.matrix(link_table_row[2]))
  OTU1abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU1),])
  OTU2abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU2),])
  distOTU1 = vegdist(cbind(OTU1abd,rep(0,length(OTU1abd))), method="bray")
  distOTU2 = vegdist(cbind(OTU2abd,rep(0,length(OTU2abd))), method="bray")

  cor1 = cor.test(distOTU1, geodist)
  cor2 = cor.test(distOTU2, geodist)

  return(c(OTU1, OTU2, cor1$estimate, cor1$p.value, cor2$estimate, cor2$p.value))
}
