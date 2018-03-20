##' This the main AIS code
##' The code will run the initial N sample of models
##' AIS is then started based on initial model weights
##' The final sample is performed using the weights that exceeded entropy threshold
##' @param filepath = parent directory above the folder with the model files 
##' @param file.name = folder where the model files and AIS will run. The final directory structure is set as paste(filepath,"/",file.name,"/",sep="")
##' @param control.name = paste(file.name,".ctl",sep="")
##' @param dat.name = paste(file.name,".dat",sep="")
##' @param tantalus = TRUE/FALSE
##' @param read.seed = TRUE/FALSE
##' @param Niter number of initial model runs
##' @param AIS.iter number of AIS model runs
##' @param final.Niter number of final model runs, recommened to be greater than N.iter and AIS.iter
##' @param m.in = c( prior mean, prior stdev) which is lognormally distributed
##' @param h.in = c( prior mean, prior stdev, lower bound, upper bound) which is a truncated beta
##' @param depl.in = c(0.50, 0.20, 0.01, 0.99,1) 
##' @param start.m.equal = FALSE
##' @param hist.yrs = the model years; seq(1916,2012,1) 
##' @param ofl.yrs = years to calculate ofl values; seq(2013,2016,1) 
##' @param depl.yr Model year associatted with the final index value, does not need to be final model year
##' @author Chantel Wetzel
##' @export

SSS.ais.fxn <- function(filepath, file.name, control.name, dat.name, 
                        tantalus=FALSE, read.seed = FALSE, entropy.level = 0.92,
                        Niter = 2000, AIS.iter = 2000, final.Niter = 5000, 
                        m.in, h.in , depl.in, start.m.equal, #hist.yrs, ofl.yrs, 
                        depl.yr)
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

 
 #Set initial seed numbers for everything-------------------------------------------------------------------------------------  
 if (read.seed == TRUE) { 
    load(paste(directory,"/seed_list",sep=""))
    seed.M      <- as.numeric(seed.list[[1]][,"seed.M"])
    seed.h      <- as.numeric(seed.list[[1]][,"seed.h"])
    seed.depl   <- as.numeric(seed.list[[1]][,"seed.depl"])
    seed.AIS    <- as.numeric(seed.list[[1]][,"seed.AIS"]) 
    seed.final  <- as.numeric(seed.list[[1]][,"seed.final"]) }
    
 if (read.seed == FALSE) {
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
   mtext("Prior Distributions", side = 3, outer=T) 
 }
 
 pdf(paste(file.name,"_priors.pdf",sep=""),width=7,height=7,)
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=",m.in[1],"sd=",m.in[3],")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",m.in[2],"sd=",m.in[4],")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",h.in[1],"sd=",h.in[2],")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",depl.in[1],"sd=",depl.in[2],")",sep=" "),main="") 
 mtext("Prior Distributions", side = 3, outer=T)
 dev.off()
 
 # Read Starter File and Change values
 starter.file <- SS_readstarter("starter.ss")
 starter.file$run_display_detail <- 0
 starter.file$detailed_age_structrure <- 2
 starter.file$parmtrace <- 0
 starter.file$cumreport <- 0
 starter.file$N_bootstraps <- 0
 starter.file$prior_like <- 0
 SS_writestarter(starter.file,file="starter.ss",overwrite=T)

 # Rename the executable if using an older executable
 if (file.exists("ss3.exe")){ file.rename("ss3.exe", "ss.exe") }
 
 print("Getting Intial Sample")
 start.time <- Sys.time()
 # Do the initial SS runs
 for (i in 1:Niter)
 {
    changeM(para=parm.vec[i,1:2])
    changeH(para=parm.vec[i,3])
    changeDepl(para=parm.vec[i,4])
        
    # Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1")}

    # Determine the model version and which files will need to be read
    # and get model dimensions
    if(i==1){
      rep.new   <- readLines("Report.sso")
      
      #Determine the SS verion
      SS_versionCode    <- rep.new[grep("#V",rep.new)]
      SS_version        <- rep.new[grep("Stock_Synthesis",rep.new)]
      SS_version        <- SS_version[substring(SS_version,1,2)!="#C"] 
      SS_versionshort   <- toupper(substr(SS_version,1,6))
      SS_versionNumeric <- as.numeric(substring(SS_versionshort,3))      
      file.name = ifelse(SS_versionNumeric >= 3.30, "ss_summary.sso", "Report.sso")

      rawrep <- read.table(file= "Report.sso" , col.names = 1:100, fill = TRUE, quote = "", 
              colClasses = "character", nrows = -1, comment.char = "")

      begin <- matchfun(string = "TIME_SERIES", obj = rawrep[,1])+2
      end   <- matchfun(string = "SPR_series",  obj = rawrep[,1])-1
      
      temptime <- rawrep[begin:end,2:3]
      endyr    <- max(as.numeric(temptime[temptime[,2]=="TIME",1]))
      startyr  <- min(as.numeric(rawrep[begin:end,2]))+2
      foreyr   <- max(as.numeric(temptime[temptime[,2]=="FORE",1]))
      
      hist.yrs <- startyr:endyr
      ofl.yrs  <- (endyr+1):foreyr
      all.yrs  <- startyr:foreyr
    }
    
    rep.new   <- readLines(file.name)
    Quant.out <- getQuant(rep.new, n=i, parm=parm.vec[i,], N=Niter, ssver = SS_versionNumeric)
    
    #Rename the report file by rep number and move to a file to save for future needs
    move.files.fxn(sim.num=i)
      
    quant.list[[1]]   <- Quant.out
    save(quant.list, file=quant.file)
 }
 end.time <- Sys.time()
 print("Intial Sample Completed Taking Total Time of:")
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
          unq.draw         <- length(unique(sir.out$get.samp))         
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
    
    effN <- 1/sum(sample.wghts*sample.wghts)
    expN <- sum(1-(1-sample.wghts)^(0.25*Niter))
    maxW <- max(sample.wghts)
    varW <- (1/length(sample.wghts))*(sum(length(sample.wghts)*sample.wghts-1)^2)
    
    entropy.out[ais]        <- entropy
    sir.list[[ais]]         <- parm.vec
    sample.wght.list[[ais]] <- sample.wghts 
        
    print(c("Entropy",entropy))
    print(c("Effective N", effN))
    print(c("Expected Unique Points", expN))
    print(c("Max Weight", maxW))
    print(c("Variance of the rescaled", varW))    
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
        if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
        if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1")}

        rep.new   <- readLines("ss_summary.sso")
        Quant.out <- getQuant(rep.new, n=i,parm = ais.parm.vec[i,], N=Niter, ssver = SS_versionNumeric)
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
    quant.good.list[[ais+1]]  <- Quant.out.good
    parm.list [[ais+1]]       <- parm.vec
    save(quant.good.list, file=quant.good.file)
    save(parm.list, file=parameters)
    
    end.time <- Sys.time()
    print("AIS Sampling Completed in Total Time:")
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
 tot.yrs <- c(hist.yrs, ofl.yrs)
 SB      <- as.data.frame(matrix(NA,nrow=length(tot.yrs),ncol=final.Niter))
 colnames(SB) = 1:final.Niter ; rownames(SB) = tot.yrs
 Bratio  <- as.data.frame(matrix(NA,nrow=(length(tot.yrs)-1),ncol=final.Niter))
 colnames(Bratio) = 1:final.Niter ; rownames(Bratio) = tot.yrs[2]:tot.yrs[length(tot.yrs)]
 TotBio  <- as.data.frame(matrix(NA,nrow=length(hist.yrs),ncol=final.Niter))
 colnames(TotBio) = 1:final.Niter ; rownames(TotBio) = hist.yrs
 OFL     <- as.data.frame(matrix(NA,nrow=length(ofl.yrs),ncol=final.Niter))
 colnames(OFL) = 1:final.Niter ; rownames(OFL) = ofl.yrs
 ForeCat <- as.data.frame(matrix(NA,nrow=length(ofl.yrs),ncol=final.Niter))
 colnames(ForeCat) = 1:final.Niter ; rownames(ForeCat) = ofl.yrs 
 
 #Do the final SS run
 for (i in 1:final.Niter)
 {
    changeM(para=final.parm.vec[i,1:2])
    changeH(para=final.parm.vec[i,3])
    changeDepl(para=final.parm.vec[i,4])
        
    #Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1") }
    
    rep.new     <- readLines(file.name)
    Quant.out   <- getQuant(rep.new, n=i,parm=final.parm.vec[i,],N=final.Niter, ssver = SS_versionNumeric)
    quant.list[[Counter+2]] <- Quant.out
    save(quant.list, file=quant.file) 
       
    RepSummary   <- RepSumFxn(rep.new, n=i, all.yrs, hist.yrs, ofl.yrs, ssver=SS_versionNumeric)
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
 rep.list[[1]]  <- TotBio[,index]
 rep.list[[2]]  <- SB[,index]
 rep.list[[3]]  <- Bratio[,index]
 rep.list[[4]]  <- OFL[,index]
 rep.list[[5]]  <- ForeCat[,index]
 
 names(rep.list) <- c("TotBio", "SB", "Bratio", "OFL", "ForeCat")
 
 parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
 colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
 
 #Get the likelihood value for each run
 likelihood    <- exp(-Quant.out.good$NLL)
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
