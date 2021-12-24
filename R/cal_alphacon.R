#' Calculating the contributions of microbial interaction on alpha diversity
#' @param env_data
#' @param  method
#' @param   ps
#' @examples
#' @return
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn}
#' @references
#' @export

cal_alphacon = function(
  env_data = data,
  method = c("forward.sel","olslm"),
  ps=ps
){

  diversity= t(estimateR(t(as.data.frame(t(vegan_otu(ps))))))

  if (method %in% c("forward.sel")) {
    res=forward.sel(diversity[,1],data,nperm=999, alpha = 0.05)
    return(res)
  }

  if (method %in% c("olslm")) {
    model=lm(diversity[,1]~.,data=as.data.frame(data))
    model.stepAIC=stepAIC(model,direction=c("both"))
    res=summary(model.stepAIC)
    return(res)
  }
}
