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
    pars.h      <- pars.truncbeta(h.in["h.start"], h.in["h.stdev"], h.in["h.LB"], h.in["h.UB"])
    alpha.h     <- pars.h[1]
    beta.h      <- pars.h[2]
    h.prior     <- dbeta(new.para$h,alpha.h,beta.h)
    
    if(depl.in["shape"]==1)
    {
      pars.depl      <- pars.truncbeta(depl.in["depl.start"], depl.in["depl.stdev"], depl.in["depl.LB"], depl.in["depl.UB"])
      alpha.depl     <- pars.depl[1]
      beta.depl      <- pars.depl[2]
      depl.prior     <- dbeta(new.para$depl,alpha.depl,beta.depl)
    }
    if(depl.in["shape"]==2){ 
      depl.prior     <- dlnorm(new.para$depl,meanlog=(log(depl.in["depl.start"])-0.5*depl.in["depl.stdev"]^2), sdlog=depl.in["depl.stdev"]) 
    }

    m.f.prior   <- dlnorm(new.para$M.f,meanlog=(log(m.in["m.f.start"])-0.5*m.in["m.f.stdev"]^2), sdlog=m.in["m.f.stdev"])
    m.m.prior   <- dlnorm(new.para$M.m,meanlog=(log(m.in["m.m.start"])-0.5*m.in["m.m.stdev"]^2), sdlog=m.in["m.m.stdev"])
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