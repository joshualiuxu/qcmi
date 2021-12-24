#' Convert the data  to phyloseq format, including OTU table and taxa table
#' @param otu_table
#' @param taxa_table
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export


trans_ps  = function(otu_table = otu_table, taxa_table =  taxa_table){

  otus=as.matrix(otu_table)
  taxa=as.matrix(taxa_table)
  taxa=taxa[row.names(otus),]
  OTU= otu_table(otus,taxa_are_rows = TRUE)
  TAX=tax_table(taxa)
  ps=phyloseq(OTU,TAX)

  return(ps)
}
