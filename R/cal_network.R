#' Correlation network calculation of microbial community data
#' @param filteredps
#' @param  N 
#' @param  r.threshold 
#' @param   p.threshold 
#' @param   lambda.min.ratio 
#' @param   nlambda 
#' @param   ncores 
#' @param   method 
#' @return 
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn} 
#' @references
#' @export


cal_network = function(ps = filteredps, N = 0,r.threshold=0.6,
                    p.threshold=0.05,lambda.min.ratio=1e-2,nlambda=20,ncores=4,method =c("spearman","pearson","spieceasi") ){
                  

if (method %in% c("spearman")) {
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
	otu_table %<>% {.[apply(., 1, sum) > 0, , drop = FALSE]}
    occor<-corAndPvalue(t(otu_table),y=NULL,use = "pairwise.complete.obs", alternative='two.sided',method='spearman')
    occor.r<-occor$cor
    occor.p<-occor$p
    occor.r[occor.p > p.threshold | abs(occor.r)<r.threshold] = 0
    ig.wgcna = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)

    
    result=list(occor.r,occor.p,method,ps_sub,ig.wgcna)
    names(result)[1] <- "occor.r"
    names(result)[2] <- "occor.p"
    names(result)[3] <- "method"
    names(result)[4] <- "phyloseq"
    names(result)[5] <- "igraph"
} 

if (method %in% c("pearson")) {
    ps_sub = filter_OTU_ps(ps = ps,Top = N)
    otu_table = as.data.frame(t(vegan_otu(ps_sub)))
	otu_table %<>% {.[apply(., 1, sum) > 0, , drop = FALSE]}
    occor<-corAndPvalue(t(otu_table),y=NULL,use = "pairwise.complete.obs", alternative='two.sided',method='pearson')
    occor.r<-occor$cor
    occor.p<-occor$p
    occor.r[occor.p > p.threshold | abs(occor.r)<r.threshold] = 0
    ig.wgcna = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)

    result=list(occor.r,occor.p,method,ps_sub,ig.wgcna)
    names(result)[1] <- "occor.r"
    names(result)[2] <- "occor.p"
    names(result)[3] <- "method"
    names(result)[4] <- "phyloseq"
    names(result)[5] <- "igraph"
} 

if (method %in% c("spieceasi")) {

    require(SpiecEasi)
    require(phyloseq)
    require(Matrix)
    ps_sub1 = filter_OTU_ps(ps = ps,Top = N)
    se.net <- spiec.easi(ps_sub1, method='mb', lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, pulsar.params=list(rep.num=50, ncores=ncores))
    ig.mb  <- adj2igraph(getRefit(se.net), vertex.attr= list(label = row.names( filteredps@otu_table)))

    sebeta <- symBeta(getOptBeta(se.net), mode='maxabs')
    elist.mb     <- summary(sebeta)


    result=list(sebeta,elist.mb ,method,ps_sub1,ig.mb,se.net)
    names(result)[1] <- "occor.r"
    names(result)[2] <- "elist"
    names(result)[3] <- "method"
    names(result)[4] <- "phyloseq"
    names(result)[5] <- "igraph"
	names(result)[6] <- "result.se"
} 

    return(result)
    
} 
