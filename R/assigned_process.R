#' Identifying ecological associations as environmental filtering and dispersal limitation
#'
#' @param link_table_row
#' @param  OTUabd
#' @param   p
#' @param   data
#' @param   cutoff
#' @param   method
#'
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

assigned_process = function(link_table_row, OTUabd, p=0.05 ,data,cutoff=0,method=c("dl", "ef")){

if (method %in% c("dl")) {

edges =link_table_row
OTUs = OTUabd
distances = data

OTUs_order = OTUs[,order(names(OTUs))]
distances_order = distances[order(row.names(distances)),]
geodist = vegdist(distances_order, method="euclid")

test_result = t(apply(edges, 1, FUN=test_link_dl, OTUabd=OTUs_order, geodist=geodist))
test_result=as.data.frame(test_result)
#test_result=cbind(row.names(test_result),test_result)
colnames(test_result)=c("OTU1","OTU2","cor1","P1","cor2","P2")
#scene1: positive links   double positive covary with geo
test_result_pick1 = test_result[which((test_result['cor1']>cutoff) & (test_result['cor2']>cutoff) & (test_result['P1']<=p) & (test_result['P2']<=p)),]

#scene2: positive links   double negative covary with geo
test_result_pick2 = test_result[which((test_result[,3]<cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]


result= rbind(test_result_pick1 ,test_result_pick2)
result_dl= result

dl=rep("yes",times=nrow(result_dl))
result_dl=cbind(result_dl,dl)

return(result_dl)
}

if (method %in% c("ef")) {

edges =link_table_row
OTUs = OTUabd
distances = data

OTUs_order = OTUs[,order(names(OTUs))]
distances_order = distances[order(row.names(distances)),]
envdist = vegdist(distances_order, method="euclid")


#scene1: positive links   double positive covary with env
edges1=edges[which((edges[,3]>=cutoff) ),]
test_result = t(apply(edges1, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
test_result=as.data.frame(test_result)
test_result_pick1 = test_result[which((test_result[,3]>cutoff) & (test_result[,5]>cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]

#scene2: positive links   double negative covary with env
edges2=edges[which((edges[,3]>=cutoff) ),]
test_result = t(apply(edges2, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
test_result=as.data.frame(test_result)
test_result_pick2 = test_result[which((test_result[,3]<cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]

#scene3: negative links   1- 1+  covary with env
edges3=edges[which((edges[,3]<cutoff) ),]
test_result = t(apply(edges3, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
test_result=as.data.frame(test_result)
test_result_pick3 = test_result[which((test_result[,3]<cutoff) & (test_result[,5] > cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]

#scene4: negative links   1+ 1- covary with env
edges4=edges[which((edges[,3]<cutoff) ),]
test_result = t(apply(edges4, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
test_result=as.data.frame(test_result)
test_result_pick4 = test_result[which((test_result[,3]>cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]


result= rbind(test_result_pick1 ,test_result_pick2,test_result_pick3,test_result_pick4)


result_ef= result

colnames(result_ef)=c("OTU1","OTU2","cor1","P1","cor2","P2")

ef=rep("yes",times=nrow(result_ef))

result_ef=cbind(result_ef,ef)

return(result_ef)
}
}

