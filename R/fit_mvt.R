##' This function calculates the covariance between the parameters and resamples based on a multivarite student's t-distribution
##' @param samp.size
##' @param para
##' @param degree
##' @author Chantel Wetzel
##' @export
fit.mvt<- function(samp.size,para,degree)
 {  
    new.dists <-list()
    
    if(start.m.equal != TRUE) {
        if (depl.in[5] == 1) { #Beta Distribution Depletion 
        para.trans    <- cbind(log(para[,1]),log(para[,2]),log((para[,3]-0.2)/(1-para[,3])), log(para[,4]/(1-para[,4])) )
        sd.log        <- apply(para.trans,2,sd)
        para.means    <- apply(para.trans,2,mean) 
        para.means    <- c(para.means[1]-0.5*sd.log[1]^2, para.means[2]-0.5*sd.log[2]^2, para.means[3], para.means[4])
        covar         <- matrix(as.numeric(cov(para.trans)),4,4)
        #Sample from the rmvt
        mult.t        <- rmvt(samp.size,sigma=covar,df=degree,method="svd")
        
        new.dists[[1]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[2]] <- exp(para.means[2] + mult.t[,2]) #M
        new.dists[[3]] <- (exp(para.means[3]+mult.t[,3])+0.20)/(1+exp(para.means[3]+mult.t[,3]))#h
        new.dists[[4]] <- exp(para.means[4]+mult.t[,4])/(1+ exp(para.means[4]+mult.t[,4])) #depletion
        }
        
        if (depl.in[5] == 2) { #Lognormal Distribution Depletion
        para.trans    <- cbind(log(para[,1]),log(para[,2]),log((para[,3]-0.2)/(1-para[,3])), log(para[,4]) )
        sd.log        <- apply(para.trans,2,sd)
        para.means    <- apply(para.trans,2,mean) 
        para.means    <- c(para.means[1]-0.5*sd.log[1]^2, para.means[2]-0.5*sd.log[2]^2, para.means[3], para.means[4]-0.5*sd.log[4]^2)
        covar         <- matrix(as.numeric(cov(para.trans)),4,4)
        #Sample from the rmvt
        mult.t        <- rmvt(samp.size,sigma=covar,df=degree,method="svd")
        
        new.dists[[1]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[2]] <- exp(para.means[2] + mult.t[,2]) #M
        new.dists[[3]] <- (exp(para.means[3]+mult.t[,3])+0.20)/(1+exp(para.means[3]+mult.t[,3]))#h
        new.dists[[4]] <- exp(para.means[4]+mult.t[,4]) #depletion
        }
    }
    
    if(start.m.equal == TRUE) {
        if (depl.in[5] == 1) { #Beta Distribution Depletion 
        para.trans    <- cbind(log(para[,1]),log((para[,2]-0.2)/(1-para[,2])), log(para[,3]/(1-para[,3])) )
        sd.log        <- apply(para.trans,2,sd)
        para.means    <- apply(para.trans,2,mean) 
        para.means    <- c(para.means[1]-0.5*sd.log[1]^2, para.means[2], para.means[3])
        covar         <- matrix(as.numeric(cov(para.trans)),3,3)
        #Sample from the rmvt
        mult.t        <- rmvt(samp.size,sigma=covar,df=degree,method="svd")
        
        new.dists[[1]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[2]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[3]] <- (exp(para.means[2]+mult.t[,2])+0.20)/(1+exp(para.means[2]+mult.t[,2])) #h
        new.dists[[4]] <- exp(para.means[3]+mult.t[,3])/(1+ exp(para.means[3]+mult.t[,3])) #depletion 
        }
        
        if (depl.in[5] == 2) { #Lognormal Distribution Depletion 
        para.trans    <- cbind(log(para[,1]),log((para[,2]-0.2)/(1-para[,2])), log(para[,3]) )
        sd.log        <- apply(para.trans,2,sd)
        para.means    <- apply(para.trans,2,mean) 
        para.means    <- c(para.means[1]-0.5*sd.log[1]^2, para.means[2], para.means[3]-0.5*sd.log[3]^2)
        covar         <- matrix(as.numeric(cov(para.trans)),3,3)
        #Sample from the rmvt
        mult.t        <- rmvt(samp.size,sigma=covar,df=degree,method="svd")
        
        new.dists[[1]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[2]] <- exp(para.means[1] + mult.t[,1]) #M
        new.dists[[3]] <- (exp(para.means[2]+mult.t[,2])+0.20)/(1+exp(para.means[2]+mult.t[,2])) #h
        new.dists[[4]] <- exp(para.means[3]+mult.t[,3]) #depletion
        }
    }

    names(new.dists) <- c("M.f","M.m", "h", "depl") 
    return(new.dists)
 }