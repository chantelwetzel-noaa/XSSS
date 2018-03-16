##' Beta Distribution for steepness
##' @param n = niter 
##' @param m = depl
##' @param s = depl.stdev
##' @param a = LB depl
##' @param b = UB depl
##' @author Chantel Wetzel
##' @export 
 rbeta.ab<- function (n, m, s, a, b)
 {   
   # CASAL's Beta
   mu    <- (m-a) / (b-a);  
   tau   <- (m-a)*(b-m)/(s^2)-1.0;
   alpha <- tau*mu;  beta <- tau*(1-mu)
   b.std <- rbeta(n, alpha, beta)
   # linear transformation from beta(0,1) to beta(a,b)
   b.out <- (b-a)*b.std + a
   
   return(b.out)
 }