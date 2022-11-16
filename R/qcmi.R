#' Quantify the strength of putative biotic associations at the community level
#' @param igraph
#' @param  OTU
#' @param   pers.cutoff
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export



qcmi = function(igraph, OTU,pers.cutoff=0){

  ag=as_adjacency_matrix(igraph,attr = 'weight')
  ag=as.matrix(ag)
  custom.cor=as.data.frame(ag)
  custom.cor.mat <- as.matrix(custom.cor)

  otu=t(OTU)

  c <- as.matrix(otu[,row.names(custom.cor)])

  rowsums.orig <- rowSums(otu)

  zero.cutoff <- ceiling(pers.cutoff * dim(c)[1]) # number of zeroes allowed in a taxon's distribution
  d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
  d <- d[rowSums(d) > 0, ]

  rel.d <- d / rowsums.orig

  obs.exp.cors.mat <- custom.cor.mat

  diag(obs.exp.cors.mat) <- 0
  # Calculate connectedness by averaging positive and negative observed - expected correlations
  connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
  connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)
  # Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
  interaction.pos <- rel.d %*% connectedness.pos
  interaction.neg <- rel.d %*% connectedness.neg

  output <- list(connectedness.neg, connectedness.pos, interaction.neg, interaction.pos)
  names(output) <- c("connectedness.neg", "connectedness.pos", "interaction.neg", "interaction.pos")

  return(output)
}
