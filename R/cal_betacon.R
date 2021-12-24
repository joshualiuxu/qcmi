#' Calculating the contributions of microbial interaction on beta diversity
#' @param env_data
#' @param  otu_abd
#' @param   method
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

cal_betacon = function(
  env_data = data,
  otu_abd=otu,
  method = c("pearson","spearman")
){

  veg.dist <- vegdist(t(otu_abd), method = 'bray')
  variable_name <- c()
  corr_res <- c()
  p_res <- c()

  for(i in 1:ncol(env_data)){
    env.dist <- vegdist(scale(env_data[, i, drop=FALSE]), "euclid")
    man1 <- mantel(veg.dist, env.dist, method =method)
    variable_name <- c(variable_name, colnames(env_data)[i])
    corr_res <- c(corr_res, man1$statistic)
    p_res <- c(p_res, man1$signif)

  }
  cor_method <- rep(method, length(p_res))
  significance <- cut(p_res, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  res_mantel <- data.frame(variable_name, cor_method, corr_res, p_res, significance)

  return(res_mantel)
}
