##' Calculate the prior probabilities for the original good parameter vectors
##' @param new.para
##' @param m.in initial M values
##' @param h.in initial steepness values
##' @param depl.in initial depletion values
##' @author Chantel Wetzel
##' @export 
 get.prior <- function(new.para, m.in, h.in, depl.in)
  {
    #Calculate the Prior(theta) values      
    # Do some calcs to help evaluate truncated beta densities later
    # bounded beta for depl, use dbeta()/cdf.up,cdf.low of fcn. to get the density value at a certain point
    pars.h      <- pars.truncbeta(h.in["h"], h.in["stdev"], h.in["LB"], h.in["UB"])
    alpha.h     <- pars.h[1]
    beta.h      <- pars.h[2]
    h.prior     <- dbeta(new.para$h,alpha.h,beta.h)
    
    if(depl.in["shape"]==1)
    {
      pars.depl      <- pars.truncbeta(depl.in["depl"], depl.in["stdev"], depl.in["LB"], depl.in["UB"])
      alpha.depl     <- pars.depl[1]
      beta.depl      <- pars.depl[2]
      depl.prior     <- dbeta(new.para$depl,alpha.depl,beta.depl)
    }
    if(depl.in["shape"]==2){ 
      depl.prior     <- dlnorm(new.para$depl,meanlog=(log(depl.in["depl"])-0.5*depl.in["stdev"]^2), sdlog=depl.in["stdev"]) 
    }

    m.f.prior   <- dlnorm(new.para$M.f,meanlog=(log(m.in["m.f"])-0.5*m.in["f.stdev"]^2), sdlog=m.in["f.stdev"])
    m.m.prior   <- dlnorm(new.para$M.m,meanlog=(log(m.in["m.m"])-0.5*m.in["m.stdev"]^2), sdlog=m.in["m.stdev"])
    if(m.in["equal.m"] == FALSE) { prior       <- m.f.prior*m.m.prior*h.prior*depl.prior }
    if(m.in["equal.m"] == TRUE)  { prior       <- m.f.prior*h.prior*depl.prior }
     
    orig.prior.list         <- list()
    orig.prior.list[[1]]    <- prior
    orig.prior.list[[2]]    <- m.f.prior
    orig.prior.list[[3]]    <- m.m.prior
    orig.prior.list[[4]]    <- h.prior
    orig.prior.list[[5]]    <- depl.prior
    names(orig.prior.list)  <- c("prior","m.f.prior","m.m.prior","h.prior","depl.prior")    

    return(orig.prior.list)
  }