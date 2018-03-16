##' Calculate the prior probabilities for the original good parameter vectors
##' @param new.para
##' @author Chantel Wetzel
##' @export 
 get.prior <- function(new.para)
  {
    #Calculate the Prior(theta) values      
    # Do some calcs to help evaluate truncated beta densities later
    # bounded beta for depl, use dbeta()/cdf.up,cdf.low of fcn. to get the density value at a certain point
    pars.h      <- pars.truncbeta(h.start, h.stdev, h.LB, h.UB)
    alpha.h     <- pars.h[1]
    beta.h      <- pars.h[2]
    h.prior     <- dbeta(new.para$h,alpha.h,beta.h)
    
    if(depl.in[5]==1)
    {
      pars.depl      <- pars.truncbeta(depl.start, depl.stdev, depl.LB, depl.UB)
      alpha.depl     <- pars.depl[1]
      beta.depl      <- pars.depl[2]
      depl.prior     <- dbeta(new.para$depl,alpha.depl,beta.depl)
    }
    if(depl.in[5]==2){ 
      depl.prior     <- dlnorm(new.para$depl,meanlog=(log(depl.start)-0.5*depl.stdev^2), sdlog=depl.stdev) 
    }

    m.f.prior   <- dlnorm(new.para$M.f,meanlog=(log(m.f.start)-0.5*m.f.stdev^2), sdlog=m.f.stdev)
    m.m.prior   <- dlnorm(new.para$M.m,meanlog=(log(m.m.start)-0.5*m.m.stdev^2), sdlog=m.m.stdev)
    if(start.m.equal != TRUE) {
      prior       <- m.f.prior*m.m.prior*h.prior*depl.prior
    }
    
    if(start.m.equal == TRUE) {
      prior       <- m.f.prior*h.prior*depl.prior
    }
     
    orig.prior.list         <- list()
    orig.prior.list[[1]]    <- prior
    orig.prior.list[[2]]    <- m.f.prior
    orig.prior.list[[3]]    <- m.m.prior
    orig.prior.list[[4]]    <- h.prior
    orig.prior.list[[5]]    <- depl.prior
    names(orig.prior.list)  <- c("prior","m.f.prior","m.m.prior","h.prior","depl.prior")    

    return(orig.prior.list)
  }