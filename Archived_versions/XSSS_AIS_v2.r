#------       XSSS        ------#
#------ Adaptive SIR V.2  ------#
#----- WO PRIOR UPDATING -------#
#------   01/02/2014      ------#

#filepath = working directory
#file.name = species names which is appended based on model runs
#m.in = c( prior mean, prior stdev) which is lognormally distributed
#h.in = c( prior mean, prior stdev, lower bound, upper bound) which is a truncated beta

#filepath="C://SSSais" 
#file.name="SHRP" 
#control.name = paste(file.name,".ctl",sep="")
#dat.name = paste(file.name,".dat",sep="")
#tantalus = TRUE/FALSE
#read.seed = TRUE/FALSE
#Niter = 1000
#AIS.iter = 1000
#final.Niter = 5000
#m.in = c(0.078, 0.077, 0.40, 0.40)
#h.in = c(0.756, 0.17, 0.20, 1)
#depl.in = c(0.50, 0.20, 0.01, 0.99,1) 
#start.m.equal = FALSE
#years = seq(1916,2012,1) 
#ofl_yrs=seq(2013,2016,1) 
#depl.yr= 2000

#Key Words for the Control File
# NatM_p_1_Fem_GP_1,  NatM_p_1_Mal_GP_1,  SR_steep
#Key Words for the Data File
# FinalDepl

SSS.ais.fxn <- function(filepath, file.name, control.name, dat.name, tantalus, read.seed, entropy.level,
                        Niter, AIS.iter, final.Niter, m.in, h.in , depl.in, start.m.equal, years, ofl_yrs, depl.yr)
{

 print("Does depl.yr match what is in the data file dummy?")

 directory  <- paste(filepath,"/",file.name,"/",sep="")
 rep.folder <- paste(directory,"report",sep="")
 save.folder<- paste(directory,"save",sep="") 
 dir.create(directory); dir.create(rep.folder); dir.create(save.folder)
 
 #Load the required packages
 require(mvtnorm)
 require(stats)
 require(r4ss)

 #Set working directory
 setwd(directory) 

 #Save the parameters
 posterior.prior.list  <- list()
 quant.list            <- list()
 quant.good.list       <- list()
 parm.list             <- list()
 rep.out.list          <- list()
 init.parms.list       <- list()
 seed.list             <- list()
 rep.list              <- list()
 
 #Output File Names
 seed.file             <- paste(save.folder,"/seed_list",sep="")
 posterior.prior.file  <- paste(save.folder,"/posterior_prior_list",sep="")
 quant.file            <- paste(save.folder,"/quant_list",sep="")
 quant.good.file       <- paste(save.folder,"/quant_good_list",sep="")
 rep.file              <- paste(save.folder,"/rep_list",sep="")
 parameters            <- paste(save.folder,"/parameters",sep="")
 init.parms.file       <- paste(save.folder,"/initial_parameters",sep="")
 entropy.file          <- paste(save.folder,"/entropy_out",sep="")
 ais.sir.file          <- paste(save.folder,"/ais_sir",sep="")
 weights.file          <- paste(save.folder,"/weights_out",sep="")

 #--------------------------------------------------------------------------------------------------------------------------
 #Beta Distribution for steepness

 # n = niter, m = depl, s = depl.stdev, a= LB depl, b = UB depl
 rbeta.ab<- function (n, m, s, a, b)
 {
   # calculate mean of corresponding standard beta dist
   #mu.std <- (m-a)/(b-a)   
   # calculate parameters of std. beta with mean=mu.std and sd=s
   #alpha <- (mu.std^2 - mu.std^3 - mu.std*s^2) / s^2
   #beta  <- (mu.std - 2*mu.std^2 + mu.std^3 - s^2 + mu.std*s^2) / s^2   
   # generate n draws from standard beta
   #b.std <- rbeta(n, alpha, beta)  
   
   # CASAL's Beta
   mu    <- (m-a) / (b-a);  
   tau   <- (m-a)*(b-m)/(s^2)-1.0;
   alpha <- tau*mu;  beta <- tau*(1-mu)
   b.std <- rbeta(n, alpha, beta)
   # linear transformation from beta(0,1) to beta(a,b)
   b.out <- (b-a)*b.std + a
   
   return(b.out)
 }
 
 #--------------------------------------------------------------------------------------------------------------------------
 # this function gets parameters of standard beta parameterization for depl and Bmsy/B0 from the truncated beta parameters
 # needed to evaluate densities in AIS
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
  
 #--------------------------------------------------------------------------------------------------------------------------
 do.sir <- function(Ncount, input, wghts)
 {
    get.samp   <- sample(1:length(wghts),size=Ncount,replace=T,prob=wghts)
    sir.vec    <- input[get.samp,]
      
    output            <- NULL
    output$sir.vec    <- sir.vec
    output$get.samp   <- get.samp
    return(output)
 } 
 
 #--------------------------------------------------------------------------------------------------------------------------

 #Calculate the prior probabilities for the original good parameter vectors
 get.prior <- function(new.para)
  {
    #Calculate the Prior(theta) values      
    # Do some calcs to help evaluate truncated beta densities later
    # bounded beta for depl, use dbeta()/cdf.up,cdf.low of fcn. to get the density value at a certain point
    pars.h      <- pars.truncbeta(h.start, h.stdev, h.LB, h.UB)
    alpha.h     <- pars.h[1]
    beta.h      <- pars.h[2]
    h.prior     <- dbeta(new.para$h,alpha.h,beta.h)
    #cdf.upper.h <- pars.h[3]
    #cdf.lower.h <- pars.h[4]
    # determine the total area under the truncated beta, this will be used to determine density later
    #dens.truncbeta.den.h <- cdf.upper.h - cdf.lower.h
    
    #cdf.upper.depl <- pars.depl[3]
    #cdf.lower.depl <- pars.depl[4]
    # determine the total area under the truncated beta, this will be used to determine density later
    #dens.truncbeta.den.depl <- cdf.upper.depl - cdf.lower.depl
    
    # evaluate the prior density for each theta (on orig. scale)
    #h.prior     <- (dbeta(new.para$h,alpha.h,beta.h)/dens.truncbeta.den.h)
    #depl.prior <- (dbeta(new.para$depl,alpha.depl,beta.depl)/dens.truncbeta.den.depl)
    
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
  
 #--------------------------------------------------------------------------------------------------------------------------  
 #Calculate the updated posterior probability for the good parameter vectors 
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
    #cdf.upper.h <- pars.h[3]
    #cdf.lower.h <- pars.h[4]
    # determine the total area under the truncated beta, this will be used to determine density later
    #dens.truncbeta.den.h <- cdf.upper.h - cdf.lower.h
    
    #cdf.upper.depl <- pars.depl[3]
    #cdf.lower.depl <- pars.depl[4]
    # determine the total area under the truncated beta, this will be used to determine density later
    #dens.truncbeta.den.depl <- cdf.upper.depl - cdf.lower.depl
    
    # evaluate the prior density for each theta (on orig. scale)
    #h.post      <- (dbeta(new.para$h, alpha.h,beta.h)/dens.truncbeta.den.h)
    #depl.post  <- (dbeta(new.para$depol, alpha.depl,beta.depl)/dens.truncbeta.den.depl)
    
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
 
 #--------------------------------------------------------------------------------------------------------------------------  
 #Calculate the new sampling weights based on the likelihood, prior, and posterior values
  get.new.wghts <- function(like,prior,posterior)
  {
    sample.wghts <- (like*prior)/posterior
    new.p        <- sample.wghts/sum(sample.wghts)
    return(new.p)
  }
 
 #------------------------------------------------------------------------------------------------------------------------------------------
 #This function calculates the covariance between the parameters and resamples based on a multivarite student's t-distribution
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
 
 #--------------------------------------------------------------------------------------------------------------------------  
 changeM <- function(para) 
 {
    para <- round(para,4)
    ctl.new <- readLines(paste(directory,control.name,sep=""))
    M.line  <- strsplit(ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)]," ")[[1]]
    M.line[c(3,4)] <- para[1]
    ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    M.line  <- strsplit(ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)]," ")[[1]]
    M.line[c(3,4)] <- para[2]
    ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    write(ctl.new,paste(directory,control.name,sep=""))
 }
 
 #--------------------------------------------------------------------------------------------------------------------------   
 changeH <- function(para)
 {
    para <- round(para,4)
    ctl.new <- readLines(paste(directory,control.name,sep=""))
    h.line  <- strsplit(ctl.new[grep("SR_steep",ctl.new)]," ")[[1]]
    h.line[c(3,4)] <- para
    ctl.new[grep("SR_steep",ctl.new)] <- paste(h.line,collapse=" ")
    write(ctl.new,paste(directory,control.name,sep=""))   
 }
 
 #--------------------------------------------------------------------------------------------------------------------------  
 changeDepl <- function(para) 
 {
    para   <- round(para,4)
    dat <- readLines(paste(directory,dat.name,sep=""))
    depl.line <- strsplit(dat[grep("FinalDepl",dat)]," ")[[1]]
    depl.line[4] <- para
    dat[grep("FinalDepl",dat)][[1]] <- paste(depl.line,collapse=" ")
    write(dat,paste(directory,dat.name,sep=""))
    
 }
 
 #--------------------------------------------------------------------------------------------------------------------------  
 #This function is currently not being used but I have kept it in here just in case I want to revert back
 getQuant <- function(n,parm,N)
 {
    if (n == 1) {
    Quant.out      <-as.data.frame(matrix(NA,nrow=N,ncol=26))
    colnames(Quant.out) <-c("M_f","M_m","h","depl","LnR0","SB0",paste("SPB_",ofl_yrs[1],sep=""),"Dep_TgtYr",
                        "Dep",paste("OFL_",ofl_yrs[2],sep=""),paste("AdjCatch_",ofl_yrs[2],sep=""),
                        "SBMSY/SB0","SSB_Unfished","SmryBio_Unfished","SSB_SPRtgt","Fstd_SPRtgt",
                        "TotYield_SPRtgt","SSB_MSY","Fstd_MSY","NLL","LL_survey","Crash","Gradient","R0_init","MissDep") }   
                                        
    rep.new         <-readLines(paste(directory,"Report.sso",sep=""))
    
    Quant.out[n,1]  <-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[n,2]  <-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[n,3]  <-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
    Quant.out[n,4]  <-parm[,4]
    Quant.out[n,5]  <-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
    Quant.out[n,6]  <-as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[n,7]  <-as.numeric(strsplit(rep.new[grep(paste("SPB_",ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[n,8]  <-as.numeric(strsplit(rep.new[grep(paste("SPB_",depl.yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[n,9]  <-as.numeric(strsplit(rep.new[grep(paste("SPB_",ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[n,10] <-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[n,11] <-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[n,12] <-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[n,13] <-as.numeric(strsplit(rep.new[grep("SSB_Unfished",rep.new)], " ")[[1]][3])
    Quant.out[n,14] <-as.numeric(strsplit(rep.new[grep("SmryBio_Unfished",rep.new)], " ")[[1]][3])
    Quant.out[n,15] <-as.numeric(strsplit(rep.new[grep("SSB_SPRtgt",rep.new)], " ")[[1]][3])
    Quant.out[n,16] <-as.numeric(strsplit(rep.new[grep("Fstd_SPRtgt",rep.new)], " ")[[1]][3])
    Quant.out[n,17] <-as.numeric(strsplit(rep.new[grep("TotYield_SPRtgt",rep.new)], " ")[[1]][3])
    Quant.out[n,18] <-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])
    Quant.out[n,19] <-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
    Quant.out[n,20] <-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    Quant.out[n,21] <-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
    Quant.out[n,22] <-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+8], " ")[[1]][2])#Crash Penalty
    Quant.out[n,23] <-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
    Quant.out[n,24] <-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
    Quant.out[n,25] <-abs(Quant.out[n,4] - Quant.out[n,8]) > 0.01
    
    return(Quant.out)
 }
 
 #Get Final Run Model Output -------------------------------------------------------------------------------------------------  
 RepSumFxn<- function(n,rep.new)
 {
    tot.yrs <- c(years,ofl_yrs)
          
    TotBio= mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste(1, x,"TIME",sep=" "),rep.new)]," ")[[1]][5]),x=years)
    SB = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste("SPB_",x,sep=""),rep.new)]," ")[[1]][3]),x=tot.yrs) 
    Bratio= mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste("Bratio_",x,sep=""),rep.new)], " ")[[1]][3]),x=tot.yrs[2]:tot.yrs[length(tot.yrs)])
    OFL = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][3]),x=ofl_yrs)
    ForeCat= mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][3]),x=ofl_yrs)    

    RepSummary <- list()
    RepSummary[[1]] <- TotBio
    RepSummary[[2]] <- SB
    RepSummary[[3]] <- Bratio
    RepSummary[[4]] <- OFL
    RepSummary[[5]] <- ForeCat

    names(RepSummary) <- c("TotBio", "SB","Bratio","OFL", "ForeCatch")
    return(RepSummary)
 }
 
 #Function that moves and saves the needed files----------------------------------------------------------------------------------------------
 move.files.fxn <- function(sim.num)
 {
    #Rename the report file by rep number and move to a file to save for future needs
    file.rename(paste(directory,"Report.sso",sep=""),paste(directory,"Report",sim.num,".sso",sep=""))
    file.copy(paste(directory,"Report",sim.num,".sso",sep=""), paste(rep.folder,"/Report",sim.num,".sso",sep=""),overwrite=T)
    file.remove(paste(directory,"Report",sim.num,".sso",sep=""))
    
    #file.copy(paste(directory,control.name,sep=""), paste(rep.folder,"/",control.name,sep=""),overwrite=T)
    #file.rename(paste(rep.folder,"/",control.name,sep=""),paste(rep.folder,"/",file.name,"_",sim.num,".ctl",sep=""))
    
    #file.copy(paste(directory,dat.name,sep=""), paste(rep.folder,"/",dat.name,sep=""),overwrite=T)
    #file.rename(paste(rep.folder,"/",dat.name,sep=""),paste(rep.folder,"/",file.name,"_",sim.num,".dat",sep=""))
 
 }
 
 #Seed Function---------------------------------------------------------------------------------------------------------------
  get.seed <- function(){ 
    Sys.sleep(1)
    t <- as.numeric(Sys.time()) ; seed <- (t-floor(t))*1e8 # set seed base upon epochs
    return(seed)
  }
 
 #Set initial seed numbers for everything-------------------------------------------------------------------------------------  
 if (read.seed == T) { 
    load(paste(directory,"/seed_list",sep=""))
    seed.M      <- as.numeric(seed.list[[1]][,"seed.M"])
    seed.h      <- as.numeric(seed.list[[1]][,"seed.h"])
    seed.depl   <- as.numeric(seed.list[[1]][,"seed.depl"])
    seed.AIS    <- as.numeric(seed.list[[1]][,"seed.AIS"]) 
    seed.final  <- as.numeric(seed.list[[1]][,"seed.final"]) }
    
 if (read.seed == F) {
    seed.M     <- get.seed()
    seed.h     <- get.seed()
    seed.depl  <- get.seed()
    seed.AIS   <- get.seed()
    seed.final <- get.seed() }
 
 #Save the seeds and write them out
 seed.list[[1]] <- cbind(seed.M,seed.h,seed.depl,seed.AIS, seed.final)
 save(seed.list, file= seed.file)
 
 #Parameter standard deviations and bounds---------------------------------------------------------------------------------
 m.f.start      <- m.in[1]      ; m.f.stdev   <- m.in[3]
 m.m.start      <- m.in[2]      ; m.m.stdev   <- m.in[4]
 h.start        <- h.in[1]      ; h.stdev     <- h.in[2]      ; h.LB     <- h.in[3]     ; h.UB     <- h.in[4]
 if(depl.in[5]==1){ depl.start    <- depl.in[1]  ; depl.stdev <- depl.in[2]  ; depl.LB <- depl.in[3] ; depl.UB <- depl.in[4]}
 if(depl.in[5]==2){ depl.start    <- depl.in[1]  ; depl.stdev <- depl.in[2] }

 #Draw M value from a lognormal distribution
 set.seed(seed.M)
 M.f <- rlnorm(Niter,meanlog=(log(m.f.start)-0.5*m.f.stdev^2), sdlog=m.f.stdev)  #should be on a log scale st dev  
 if(start.m.equal == TRUE) { set.seed(seed.M) }
 M.m <- rlnorm(Niter,meanlog=(log(m.m.start)-0.5*m.m.stdev^2), sdlog=m.m.stdev)

 #Draw h from beta distribution
 set.seed(seed.h)
 h <- rbeta.ab(Niter,h.start, h.stdev, h.LB, h.UB)    
 
 #Draw depl from beta distribution
 set.seed(seed.depl)
 if(depl.in[5]==1) {depl <- rbeta.ab(Niter,depl.start, depl.stdev, depl.LB, depl.UB)}
 if(depl.in[5]==2) {depl <- rlnorm(Niter,meanlog=(log(depl.start)-0.5*depl.stdev^2), sdlog=depl.stdev)}
   
 parm.vec <- cbind.data.frame(M.f, M.m, h, depl)
 parm.vec <- round(parm.vec,4)
 
 if(tantalus == FALSE) {
   windows(record=TRUE)
   par(mfrow=c(2,2),oma=c(3,3,4,3))
   hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=",m.in[1],"sd=",m.in[3],")", sep=" "), main="")
   hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",m.in[2],"sd=",m.in[4],")", sep=" "), main="")
   hist(parm.vec[,3],xlab= paste("Steepness (mean=",h.in[1],"sd=",h.in[2],")", sep=" "),main="")
   hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr,"(mean=",depl.in[1],"sd=",depl.in[2],")",sep=" "),main="") 
   mtext("Prior Distributions", side = 3, outer=T) }
 
 pdf(paste(file.name,"_priors.pdf",sep=""),width=7,height=7,)
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=",m.in[1],"sd=",m.in[3],")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",m.in[2],"sd=",m.in[4],")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",h.in[1],"sd=",h.in[2],")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",depl.in[1],"sd=",depl.in[2],")",sep=" "),main="") 
 mtext("Prior Distributions", side = 3, outer=T)
 dev.off()
 
 #Read Starter File and Change values
 starter.file <- SS_readstarter(paste(directory,"starter.ss",sep=""))
 starter.file$run_display_detail <- 0
 starter.file$detailed_age_structrure <- 0
 starter.file$parmtrace <- 0
 starter.file$cumreport <- 0
 starter.file$N_bootstraps <- 0
 starter.file$prior_like <- 0
 SS_writestarter(starter.file,file="starter.ss",overwrite=T)
 
 print("Getting Intial Sample")
 start.time <- Sys.time()
 #Do the initial SS runs
 for (i in 1:Niter)
 {
    changeM(para=parm.vec[i,1:2])
    changeH(para=parm.vec[i,3])
    changeDepl(para=parm.vec[i,4])
        
    #Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS3 -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss3.exe -nohess > out.txt 2>&1")}
    
    Quant.out <- getQuant(n=i,parm=parm.vec[i,],N=Niter)
    
    #Rename the report file by rep number and move to a file to save for future needs
    move.files.fxn(sim.num=i)
      
    quant.list[[1]]   <- Quant.out
    save(quant.list, file=quant.file)
 }
 end.time <- Sys.time()
 print(end.time - start.time)
 
 
 #Remove runs where depletion was not met
 Quant.out.good <- Quant.out[Quant.out$MissDep < 1 & is.na(Quant.out$MissDep)!=T & Quant.out$Crash == 0,]
  
 init.parms.list[[1]]   <- parm.vec #Save the initial parameter before removing them
 parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
 colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
 
 #Get the likelihood value for each run
 #perhaps change this to exp(-LLsurvey) instead or exp(-NLL) for total likelihood value
 likelihood <- exp(-Quant.out.good$NLL)
 
 #Calculate the prior probability based on the initial good parameter vectors
 prior.all     <- get.prior(new.para=parm.vec) 

 #Save the parameters 
 posterior.prior.list[[1]]    <- prior.all
 quant.good.list[[1]]         <- Quant.out.good
 parm.list[[1]]               <- parm.vec   
 
 #Write saved files with year index
 save(quant.good.list, file=quant.good.file)
 save(posterior.prior.list, file=posterior.prior.file)
 save(parm.list, file=parameters)
 save(init.parms.list, file = init.parms.file)
 
 #Set seed for multiple sampling
 set.seed(seed.AIS)
 entropy.out <- rep(1,10)
 sir.list <- list()
 sample.wght.list <- list()
 
 Niter = AIS.iter
 #Run AIS up to 10x or when entropy >= 0.92
 for (ais in 1:10)
 {  
    start.time <- Sys.time()
    Counter = ais       
    #The first sample weights are proportional to the likelihood
    if (ais == 1)
      {        
          prior            <- prior.all$prior
          posterior.all    <- get.prior(new.para=parm.vec)
          p                <- likelihood/sum(likelihood) 
          sample.wghts     <- p 
          sir.out          <- do.sir(Ncount=0.25*Niter,input= parm.vec, wghts=sample.wghts)
          #Check for how many unique SIR draws 
          unq.draw           <- length(unique(sir.out$get.samp))         
      }
    
    #Here after sample weights are = likelihood*prior/pr
    if (ais > 1 )
      {
         ais.parm.vec       <- parm.vec
         p                  <- likelihood/sum(likelihood)
         posterior.all      <- get.new.posterior(new.para=ais.parm.vec,degree=5)
         posterior          <- posterior.all$posterior
         prior.all          <- get.prior(new.para=ais.parm.vec) 
         prior              <- prior.all$prior
         sample.wghts       <- get.new.wghts(p, prior, posterior)
         sir.out            <- do.sir(Ncount=0.25*Niter,input= ais.parm.vec, wghts=sample.wghts)       
         #Check for how many unique SIR draws 
         unq.draw           <- length(unique(sir.out$get.samp))
      }
    
    #Store the SIR trajectories in a data frame and name columns
    sir.vec        <- as.data.frame(sir.out$sir.vec)
    names(sir.vec) <- c("M.f","M.m", "h", "depl")

    posterior.prior.list[[ais+1]]  <- c(prior.all, posterior.all) 
    
    #Calculte the entropy West 1993
    if (all(sample.wghts>0)) { entropy <- -sum(sample.wghts *(log(sample.wghts)/log(length(sample.wghts)))) 
    }  else { entropy <- -sum(sample.wghts[sample.wghts>0] * (log(sample.wghts[sample.wghts>0])/log(length(sample.wghts)))) }
    
    entropy.out[ais]        <- entropy
    sir.list[[ais]]         <- parm.vec
    sample.wght.list[[ais]] <- sample.wghts 
        
    print(entropy)    
    if (entropy >=entropy.level) break()
                         
    #This is where sample from the new parameter values
    new.dists           <- fit.mvt(Niter,para=sir.vec,degree=5)
    ais.parm.vec        <- cbind.data.frame(new.dists$M.f,new.dists$M.m,new.dists$h,new.dists$depl)
    names(ais.parm.vec) <- c("M.f","M.m", "h","depl")
    
    #Do the initial SS runs
    for (i in 1:Niter)
    {
        changeM(para=ais.parm.vec[i,1:2])
        changeH(para=ais.parm.vec[i,3])
        changeDepl(para=ais.parm.vec[i,4])
        
        #Run Simple Stock Synthesis
        if (tantalus == TRUE) { system("./SS3 -nohess > out.txt 2>&1")  }
        if (tantalus == FALSE){ shell("ss3.exe -nohess > out.txt 2>&1")}

        Quant.out <- getQuant(n=i,parm = ais.parm.vec[i,],N=Niter)
        quant.list[[ais+1]] <- Quant.out
        save(quant.list, file=quant.file)  
    }
    
    #Remove runs Crashed Runs and where Depletion was not Met
    Quant.out.good <- Quant.out[Quant.out$MissDep < 1 & is.na(Quant.out$MissDep)!=T & Quant.out$Crash == 0,]
 
    parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
    colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
    
    #Get the likelihood value for each run
    #trans.like <- exp(-(Quant.out.good$LL_survey-min(Quant.out.good$LL_survey)))
    likelihood <- exp(-Quant.out.good$NLL)

    #Save the quantities from the report files
    quant.good.list[[ais+1]]   <- Quant.out.good
    parm.list [[ais+1]]   <- parm.vec
    save(quant.good.list, file=quant.good.file)
    save(parm.list, file=parameters)
    
    end.time <- Sys.time()
    print(end.time - start.time)
    
 }#End AIS loop
 
 if (entropy.out[10] < 0.92) {
    xx <- sort(entropy.out,index.return=T)$ix[10]
    parm.vec <- sir.list[[xx]]
    sample.wghts <- sample.wght.list[[xx]]
    
    save(entropy.out, file = entropy.file)
    save(sir.list, file= ais.sir.file)
    save(sample.wght.list, file=weights.file)
 }  

 set.seed(seed.final)
 print("Getting Final Sample")
 start.time <- Sys.time()
 #Create the final parameter distributions
 final.sir        <- do.sir(Ncount=final.Niter, input=parm.vec, wghts=sample.wghts)
 final.sir.vec    <- as.data.frame(final.sir$sir.vec)
 names(final.sir.vec)<- c("M.f","M.m", "h", "depl")

 final.parm.vec        <- cbind.data.frame(final.sir.vec$M.f,final.sir.vec$M.m,final.sir.vec$h,final.sir.vec$depl)
 names(final.parm.vec) <- c("M.f","M.m", "h", "depl")
 
 #Create Storage matrixes
 tot.yrs <- c(years,ofl_yrs)
 SB      <- as.data.frame(matrix(NA,nrow=length(tot.yrs),ncol=final.Niter))
 colnames(SB) = 1:final.Niter ; rownames(SB) = tot.yrs
 Bratio  <- as.data.frame(matrix(NA,nrow=(length(tot.yrs)-1),ncol=final.Niter))
 colnames(Bratio) = 1:final.Niter ; rownames(Bratio) = tot.yrs[2]:tot.yrs[length(tot.yrs)]
 TotBio  <- as.data.frame(matrix(NA,nrow=length(years),ncol=final.Niter))
 colnames(TotBio) = 1:final.Niter ; rownames(TotBio) = years
 OFL     <- as.data.frame(matrix(NA,nrow=length(ofl_yrs),ncol=final.Niter))
 colnames(OFL) = 1:final.Niter ; rownames(OFL) = ofl_yrs
 ForeCat <- as.data.frame(matrix(NA,nrow=length(ofl_yrs),ncol=final.Niter))
 colnames(ForeCat) = 1:final.Niter ; rownames(ForeCat) = ofl_yrs 
 
 #Do the final SS run
 for (i in 1:final.Niter)
 {
    changeM(para=final.parm.vec[i,1:2])
    changeH(para=final.parm.vec[i,3])
    changeDepl(para=final.parm.vec[i,4])
        
    #Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS3 -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss3.exe -nohess > out.txt 2>&1") }
    rep.new      <- readLines(paste(directory,"Report.sso",sep=""))
    
    Quant.out <- getQuant(n=i,parm=final.parm.vec[i,],N=final.Niter)
    quant.list[[Counter+2]] <- Quant.out
    save(quant.list, file=quant.file) 
       
    RepSummary <-  RepSumFxn(n=i,rep.new)
    SB[,i]       <- RepSummary$SB
    TotBio[,i]   <- RepSummary$TotBio
    Bratio[,i]   <- RepSummary$Bratio
    OFL[,i]      <- RepSummary$OFL
    ForeCat[,i]  <- RepSummary$ForeCat
     
    #Rename the report file by rep number and move to a file to save for future needs
    move.files.fxn(sim.num=i) 
    
 }
 end.time <- Sys.time()
 print(end.time - start.time)
 
 #Remove runs Crashed Runs and where Depletion was not Met
 index  <- (Quant.out$MissDep == FALSE & is.na(Quant.out$MissDep)!=T & Quant.out$Crash == 0)
 
 Quant.out.good <- Quant.out[index,]
 rep.list[[1]] <- TotBio[,index]
 rep.list[[2]] <- SB[,index]
 rep.list[[3]] <- Bratio[,index]
 rep.list[[4]] <- OFL[,index]
 rep.list[[5]] <- ForeCat[,index]
 
 parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
 colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
 
 #Get the likelihood value for each run
 likelihood <- exp(-Quant.out.good$NLL)
 prior.all     <- get.prior(new.para=parm.vec[,1:4])  
 posterior.all <- get.new.posterior(new.para=parm.vec[,1:4],degree=5)
 
 #Save the final trajectories, liklihood values, and probabilities
 posterior.prior.list[[Counter+2]]  <- c(prior.all, posterior.all)
 quant.good.list[[Counter+2]]   <- Quant.out.good
 parm.list[[Counter+2]] <- parm.vec

 #Write saved files with year index
 save(quant.good.list, file=quant.good.file)
 save(rep.list, file= rep.file)
 save(posterior.prior.list, file=posterior.prior.file)
 save(parm.list, file=parameters)
 
 if (tantalus == FALSE) {
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=", round(mean(parm.vec[,1]),2),"sd=",round(sd(parm.vec[,1]),2),")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",round(mean(parm.vec[,2]),2),"sd=",round(sd(parm.vec[,2]),2),")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",round(mean(parm.vec[,3]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",round(mean(parm.vec[,4]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 mtext("Posterior Distributions", side = 3, outer=T) }
 
 pdf(paste(file.name,"_posteriors.pdf",sep=""),width=7,height=7,)
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=", round(mean(parm.vec[,1]),2),"sd=",round(sd(parm.vec[,1]),2),")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",round(mean(parm.vec[,2]),2),"sd=",round(sd(parm.vec[,2]),2),")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",round(mean(parm.vec[,3]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",round(mean(parm.vec[,4]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 mtext("Posterior Distributions", side = 3, outer=T)
 dev.off()
 
 return("Boom Pow Assessed!")
}
