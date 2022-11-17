#' Random Matrix Theory to filter correlation coefficient
#' @param cormat
#' @param  lcor 
#' @param  hcor 
#' @return 
#' @author Contact: Xu Liu \email{xliu@@issas.ac.cn} 
#' @references
#' @export

#RMT
#following by an RMT-based approach that determines the correlation cut-off threshold in an automatic fashion. 
#Random matrix theory was initially proposed in the 1960s as a procedure to identify phase transitions associated with noise in physics and material science,
#and was later adopted for studying the behaviours of many other complex systems, including gene co-expression network construction for
#predicting gene functions and molecular ecological network construction.	
rmt = function(cormat,lcor=0.6, hcor=0.8){
    s <- seq(0, 3, 0.1)
    pois <- exp(-s)
    geo <- 0.5 * pi * s * exp(-0.25 * pi * s^2)
    ps <- NULL  
    for(i in seq(lcor, hcor, 0.01)){
        cormat1 <- abs(cormat)
        cormat1[cormat1 < i] <- 0  
        eigen <- sort(eigen(cormat1)$value)
        ssp <- smooth.spline(eigen, control.spar = list(low = 0, high = 3)) 
        nnsd1 <- density(nnsd(ssp$y))
        nnsdpois <- density(nnsd(pois))
        chival1 <- sum((nnsd1$y - nnsdpois$y)^2/nnsdpois$y/512)
        ps <- rbind(ps, chival1)
        print(i*100)
    }
    ps <- cbind(ps, c(seq(lcor, hcor, 0.01)))
    tc <- ps[ps[,1] == min(ps[,1]), 2]
    tc
}