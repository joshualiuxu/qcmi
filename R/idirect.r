#' Disentangling direct from indirect relationships in association networks
#' @param pathname
#' @param  filename 
#' @param  outputname 
#' @return 
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn} 
#' @references
#' @export

#idirect
#idirect addresses the three main problems of ill-conditioning, self-loop and interaction strength overflow that may exist in previous methods 

idirect = function(pathname=NULL, filename = NULL, outputname = NULL){

os <- import("os")

os$chdir(pathname)

time =import("time") 

idir =import("idirect") 

nh =import("net_handler") 

fh =import("file_handler") 

t0 = time$time()

G = fh$read_file_weighted_edges(filename, t0)

S = idir$direct_association(G[[1]], t0=t0)

S2 = nh$merge(G[[1]], S[[1]])

St = fh$save_sorted_turple(S2, in_file=filename)

fh$save_file_weighted_edges(St, outputname, t0=t0)

}