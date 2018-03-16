##' Calculate the updated posterior probability for the good parameter vectors 
##' @param new.para
##' @param degree
##' @author Chantel Wetzel
##' @export 
 get.new.posterior <- function(new.para,degree)
  {
    #Calculate the new means and standard deviations 
    mean.h    <- mean(new.para$h)      ; h.sd         <- sd(new.para$h)
    mean.depl <- mean(new.para$depl)   ; depl.sd      <- sd(new.para$depl)
    mean.m.f  <- mean(new.para$M.f)    ; m.f.sd       <- sd(new.para$M.f)
    mean.m.m  <- mean(new.para$M.f)    ; m.m.sd       <- sd(new.para$M.m)
    
    #Calculate the new Pr(theta) values that are NOT used in the denominator but for final comparison of dists
    pars.h      <- pars.truncbeta(mean.h, h.sd, h.LB, h.UB)
    alpha.h     <- pars.h[1]
    beta.h      <- pars.h[2]
    h.post      <- dbeta(new.para$h,alpha.h,beta.h)
    
    if(depl.in[5]==1)
    {
      pars.depl      <- pars.truncbeta(mean.depl, depl.sd, depl.LB, depl.UB)
      alpha.depl     <- pars.depl[1]
      beta.depl      <- pars.depl[2]
      depl.post      <- dbeta(new.para$depl, alpha.depl, beta.depl)
    }
    
    if(depl.in[5]==2){
      depl.post      <- dlnorm(new.para$depl,meanlog=(log(mean.depl)-0.5*depl.sd^2), sdlog=depl.sd)
    }

    m.f.post    <- dlnorm(new.para$M.f,meanlog=(log(mean.m.f)-0.5*m.f.sd^2), sdlog=m.f.sd)
    m.m.post    <- dlnorm(new.para$M.m,meanlog=(log(mean.m.m)-0.5*m.m.sd^2), sdlog=m.m.sd)
    
    #Calculate the Pr(theta) values based on the student t distribution
    if (start.m.equal != TRUE) {
      covar.x     <- matrix(cov(new.para),4,4)
      mean.vec    <- c(mean.m.f, mean.m.m, mean.h,mean.depl)
      posterior   <- numeric(dim(new.para)[1])
        for(i in 1:(dim(new.para)[1])) { 
          posterior[i]   <- dmvt(new.para[i,1:4],delta=mean.vec,sigma=covar.x,df=degree,type='shifted',log=FALSE) }
    }
    
    #Calculate the Pr(theta) values based on the student t distribution
    if (start.m.equal == TRUE) {
      covar.x     <- matrix(cov(new.para[,c(1,3,4)]),3,3)
      mean.vec    <- c(mean.m.f, mean.h, mean.depl)
      posterior   <- numeric(dim(new.para)[1])
        for(i in 1:(dim(new.para)[1])) { 
          posterior[i]   <- dmvt(new.para[i,c(1,3,4)],delta=mean.vec,sigma=covar.x,df=degree,type='shifted',log=FALSE) }
    }
    
    new.post.list         <- list()
    new.post.list[[1]]    <- posterior
    new.post.list[[2]]    <- m.f.post
    new.post.list[[3]]    <- m.m.post
    new.post.list[[4]]    <- h.post
    new.post.list[[5]]    <- depl.post
   
    names(new.post.list)  <- c("posterior","m.f.post","m.m.post","h.post","depl.post")
    return(new.post.list)  
  } 