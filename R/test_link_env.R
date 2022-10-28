#' Classifying ecological associations of environmental filtering
#' @param   link_table_row
#' @param   OTUabd
#' @param   envdist
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export


test_link_env = function(link_table_row, OTUabd, envdist){
	OTU1 = as.character(as.matrix(link_table_row[1]))
	OTU2	 = as.character(as.matrix(link_table_row[2]))
	OTU1abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU1),])
	OTU2abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU2),])
	distOTU1 = vegdist(cbind(OTU1abd,rep(0,length(OTU1abd))), method="bray")
	distOTU2 = vegdist(cbind(OTU2abd,rep(0,length(OTU2abd))), method="bray")

	cor1 = cor.test(distOTU1, envdist)
	cor2 = cor.test(distOTU2, envdist)

	return(c(OTU1, OTU2, cor1$estimate, cor1$p.value, cor2$estimate, cor2$p.value))
}
