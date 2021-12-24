#' Correlation network calculation of microbial community data
#' @param filteredps
#' @param  N
#' @param  r.threshold
#' @param   p.threshold
#' @param   method
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export


cal_network = function(ps = filteredps, N = 0,r.threshold=0.6,
                       p.threshold=0.05,method =c("spearman","pearson") ){


  if (method %in% c("spearman")) {
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
    occor<-corAndPvalue(t(otu_table),y=NULL,use = "pairwise.complete.obs", alternative='two.sided',method='spearman')
    occor.r<-occor$cor
    occor.p<-occor$p
  }

  if (method %in% c("pearson")) {
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
    occor<-corAndPvalue(t(otu_table),y=NULL,use = "pairwise.complete.obs", alternative='two.sided',method='pearson')
    occor.r<-occor$cor
    occor.p<-occor$p
  }

  occor.r[occor.p > p.threshold | abs(occor.r)<r.threshold] = 0
  ig.wgcna = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)

  result=list(occor.r,occor.p,method,ps_sub,ig.wgcna)
  names(result)[1] <- "occor.r"
  names(result)[2] <- "occor.p"
  names(result)[3] <- "method"
  names(result)[4] <- "phyloseq"
  names(result)[5] <- "igraph"

  return(result)

}
