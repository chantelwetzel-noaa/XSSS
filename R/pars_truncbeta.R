##' this function gets parameters of standard beta parameterization for depl and Bmsy/B0 from the truncated beta parameters
##' needed to evaluate densities in AIS
##' @param m  
##' @param s
##' @param low
##' @param up
##' @author Chantel Wetzel
##' @export 
 pars.truncbeta <- function(m,s,low,up)
  {  
    mu  <- (m-low)/(up-low)
    # calculate parameters of std. beta with mean=mu.std and s
    #alpha       <- as.numeric((mu.std^2 - mu.std^3 - mu.std*s^2) / s^2)
    #beta        <- as.numeric((mu.std - 2*mu.std^2 + mu.std^3 - s^2 + mu.std*s^2) / s^2)
    tau <- (m-low)*(up-m)/(s^2)-1.0;
    alpha <- tau*mu;  beta <- tau*(1-mu)
    cdf.up      <- pbeta(as.numeric(up),alpha,beta)
    cdf.low     <- pbeta(as.numeric(low),alpha,beta)
    return(c(alpha,beta,cdf.up,cdf.low))
  }