##' This the main AIS code
##' The code will run the initial N sample of models
##' AIS is then started based on initial model weights
##' The final sample is performed using the weights that exceeded entropy threshold
##' @param filepath = parent directory above the folder with the model files 
##' @param control.name 
##' @param dat.name
##' @param tantalus = TRUE/FALSE,  Option for use with linux systems, NWFSC specific
##' @param read.seed = TRUE/FALSE, Option to allow the user to predefine seeds from previous model run, reads the seed_list file
##' @param Niter number of initial model runs
##' @param AIS.iter number of AIS model runs
##' @param final.Niter number of final model runs, recommened to be greater than N.iter and AIS.iter
##' @param m.in = c(prior mean f, prior mean m, prior stdev, prior stdev, start m equal) lognormally distributed
##' @param h.in = c(prior mean, prior stdev, lower bound, upper bound, distribution (trunacted beta only))
##' @param depl.in = c(0.50, 0.20, 0.01, 0.99, distribution (1 = truncated beta, 2 = lognormal)) 
##' @param fmsy.m.in NOT YET FULLY IMPLEMENTED
##' @param bmsy.b0.in NOT YET FULLY IMPLEMENTED
##' @author Chantel Wetzel
##' @export
##' @seealso \code{\link{rbeta_ab_fn}}, \code{\link{pars_truncbeta_fn}},
##' \code{\link{do_sir_fn}}, \code{\link{get_prior_fn}},
##' \code{\link{get_new_posteriors_fn}}, \code{\link{get_new_wghts_fn}},
##' \code{\link{fit_mvt_fn}}, \code{\link{change_m_fn}},
##' \code{\link{change_h_fn}}, \code{\link{change_depl_fn}},
##' \code{\link{get_quant_fn}}, \code{\link{summary_fn}},
##' \code{\link{move_file_fn}}, \code{\link{get_see_fn}},
##' \code{\link{match_fn}}, \code{\link{define_matrix_fn}}
##' \cod{\linkd{plot_fn}}
##' @import mvtnorm
##' @import stats
##' @import r4ss

SSS.ais.fxn <- function(filepath, control.name, dat.name, 
                        tantalus=FALSE, read.seed = FALSE, entropy.level = 0.92,
                        Niter = 2000, AIS.iter = 2000, final.Niter = 5000, 
                        m.in, h.in , depl.in, fmsy.m.in = NULL, bmsy.b0.in = NULL) 
                        
{

 #Load the required packages
 require(mvtnorm)
 require(stats)
 require(r4ss)

 directory  <- paste(filepath,"/run/",sep="")
 rep.folder <- paste(directory,"report",sep="")
 save.folder<- paste(directory,"save",sep="") 
 dir.create(directory, showWarnings = FALSE)
 dir.create(rep.folder, showWarnings = FALSE)
 dir.create(save.folder, showWarnings = FALSE)
 
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
    seed.M       <- as.numeric(seed.list[[1]][,"seed.M"])
    seed.h       <- as.numeric(seed.list[[1]][,"seed.h"])
    seed.depl    <- as.numeric(seed.list[[1]][,"seed.depl"])
    seed.fmsy.m  <- as.numeric(seed.list[[1]][,"seed.fmsy.m"])
    seed.bmsy.b0 <- as.numeric(seed.list[[1]][,"seed.bmsy.b0"])
    seed.AIS     <- as.numeric(seed.list[[1]][,"seed.AIS"]) 
    seed.final   <- as.numeric(seed.list[[1]][,"seed.final"]) }
    
 if (read.seed == FALSE) {
    seed.M       <- get.seed()
    seed.h       <- get.seed()
    seed.depl    <- get.seed()
    seed.fmsy.m  <- get.seed()
    seed.bmsy.b0 <- get.seed()
    seed.AIS     <- get.seed()
    seed.final   <- get.seed() }
 
 #Save the seeds and write them out
 seed.list[[1]] <- cbind(seed.M,seed.h,seed.depl, seed.fmsy.m, seed.bmsy.b0, seed.AIS, seed.final)
 save(seed.list, file= seed.file)
 
 #Parameter standard deviations and bounds---------------------------------------------------------------------------------
 names(m.in)    = c("m.f", "m.m", "f.stdev", "m.stdev", "equal.m")
 names(h.in)    = c("h", "stdev", "LB", "UB")
 names(depl.in) = c("depl", "stdev", "LB", "UB", "shape") 

 #Draw M value from a lognormal distribution
 set.seed(seed.M)
 M.f <- rlnorm(Niter,meanlog=(log(m.in["m.f"])-0.5*m.in["f.stdev"]^2), sdlog=m.in["f.stdev"])  #should be on a log scale st dev  
 if(m.in["equal.m"] == TRUE) { set.seed(seed.M) }
 M.m <- rlnorm(Niter,meanlog=(log(m.in["m.m"])-0.5*m.in["m.stdev"]^2), sdlog=m.in["m.stdev"])

 #Draw h from beta distribution
 set.seed(seed.h)
 h <- rbeta.ab(Niter,h.in["h"], h.in["stdev"], h.in["LB"], h.in["UB"])    
 
 #Draw depl from beta or lognormal distribution
 set.seed(seed.depl)
 if(depl.in["shape"]==1) {depl <- rbeta.ab(Niter,depl.in["depl"], depl.in["stdev"], depl.in["LB"], depl.in["UB"])}
 if(depl.in["shape"]==2) {depl <- rlnorm(Niter,meanlog=(log(depl.in["depl"])-0.5*depl.in["stdev"]^2), sdlog=depl.in["stdev"])}

 #Draw fmsy/m from beta or lognormal distribution
 if (!is.null(fmsy.m.in)){
  names(fmsy.m.in) = c("fmsy.m", "stdev", "LB", "UB", "shape")
  set.seed(seed.fmsy.m)
  if(fmsy.m.in["shape"]==1) {fmsy.m <- rbeta.ab(Niter,fmsy.m.in["fmsy.m"], fmsy.m.in["stdev"], fmsy.m.in["LB"], fmsy.m.in["UB"])}
  if(fmsy.m.in["shape"]==2) {fmsy.m <- rlnorm(Niter,meanlog=(log(fmsy.m.in["fmsy.m"])-0.5*fmsy.m.in["stdev"]^2), sdlog=fmsy.m.in["stdev"])}
 }

 if (is.null(fmsy.m.in) & !is.null(bmsy.b0.in)) {c("Need to define bmsy.b0 since fmsy.m is being used")}
 if (!is.null(bmsy.b0.in) ) {
  names(fmsy.m.in) = c("bmsy.b0", "stdev", "LB", "UB", "shape")
  set.seed(seed.bmsy.b0)
  if(bmsy.b0.in["shape"]==1) {bmsy.b0 <- rbeta.ab(Niter,fmsy.m.in["bmsy.b0"], fmsy.m.in["stdev"], fmsy.m.in["LB"], fmsy.m.in["UB"])}
  if(bmsy.b0.in["shape"]==2) {bmsy.b0 <- rlnorm(Niter,meanlog=(log(fmsy.m.in["bmsy.b0"])-0.5*fmsy.m.in["stdev"]^2), sdlog=fmsy.m.in["stdev"])}
 }


 parm.vec <- cbind.data.frame(M.f, M.m, h, depl)
 parm.vec <- round(parm.vec,4)

 # Move model files into the species folder
 file.copy(paste0(filepath,"/", "starter.ss"), "starter.ss")
 file.copy(paste0(filepath,"/", "forecast.ss"), "forecast.ss")
 file.copy(paste0(filepath,"/", control.name), control.name)
 file.copy(paste0(filepath,"/", dat.name), dat.name)
 file.copy(paste0(filepath,"/", "ss3.exe"), "ss.exe")
 file.copy(paste0(filepath,"/", "ss.exe"),  "ss.exe")

 # Find depletion year
 dat     <- readLines(dat.name)
 depl.yr <- as.numeric(strsplit(dat[grep("FinalDepl",dat)],"[[:blank:]]+")[[1]][1])

 # Create Prior Distribution Plots
 pdf(paste0(save.folder,"/priors.pdf"),width=7,height=7,)
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=",m.in[1],"sd=",m.in[3],")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",m.in[2],"sd=",m.in[4],")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",h.in[1],"sd=",h.in[2],")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",depl.in[1],"sd=",depl.in[2],")",sep=" "),main="") 
 mtext("Prior Distributions", side = 3, outer=T)
 dev.off()

 file.name = "ss_summary.sso"
 
 # Read Starter File and Change values
 starter.file <- SS_readstarter("starter.ss", verbose = FALSE)
 starter.file$datafile <- dat.name
 starter.file$ctlfile <- control.name
 starter.file$run_display_detail <- 0
 starter.file$detailed_age_structure <- 2
 starter.file$last_estimation_phase <- 5
 starter.file$parmtrace <- 0
 starter.file$cumreport <- 0
 starter.file$N_bootstraps <- 0
 starter.file$prior_like <- 0
 starter.file$jitter <- 0 
 starter.file$SPR_basis <- 4 # Report 1-SPR 
 SS_writestarter(starter.file,file="starter.ss",overwrite=T, verbose=FALSE, warn=FALSE)

 data.file <- SS_readdat_3.30(dat.name, verbose = FALSE)
 data.file$CPUEinfo[,"SD_Report"] = 1
 names = c(rep("#", dim(data.file$CPUE)[1] -1), "#FinalDepl")
 data.file$CPUE = cbind(data.file$CPUE, names)
 SS_writedat_3.30(data.file, outfile = dat.name, overwrite = T, verbose = FALSE)
 ind = data.file$fleetinfo[,"type"]
 ind = ind == 3
 survey.names = data.file$fleetnames[ind]
 temp = 1:length(ind)
 survey.fleet.number = temp[ind]
 survey.yrs = NULL
 for (a in 1:length(survey.names)){
  find = data.file$CPUE$index == survey.fleet.number[a]
  survey.yrs = c(survey.yrs, as.numeric(data.file$CPUE$year[find]))
 }

 # Run Simple Stock Synthesis
 if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
 if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1")}

 # Determine the model version and which files will need to be read
 # and get model dimensions
 rep.new   <- readLines("Report.sso")
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

 # Determine the number of fleets and the catch
 begin <- matchfun(string = "CATCH", obj = rawrep[,1])+2
 end   <- matchfun(string = "TIME_SERIES", obj = rawrep[,1])-1
 temp  <- rawrep[begin:end, 1:18]
 names <- unique(temp[,2]) # This is a list of fishery names with catch associated with them
 fleet.num <- unique(temp[,1])

 total.dead = total.catch = 0 
 for (i in 1:length(fleet.num)){
   killed = mapply(function(x) killed = as.numeric(strsplit(rep.new[grep(paste(fleet.num[i], names[i], x, sep=" "),rep.new)]," ")[[1]][14]), x = hist.yrs)
   total.dead = total.dead + killed
 }
  
 Catch = total.dead
 names(Catch) = hist.yrs

 # Determine how many suvey fleets are included
 find = matchfun(string = "Surv_like", obj = rawrep[,1])
 temp = as.numeric(rawrep[find,3:10])
 survey.list = temp[!is.na(temp)]
 n.survey = sum(survey.list != 0 ) - 1 # Substract 1 to remove depl survey


 # Determine if the model has added variance is used in the control file
 rawctl <- read.table(file= control.name , col.names = 1:20, fill = TRUE, quote = "", 
        colClasses = "character", nrows = -1, comment.char = "")

 for(a in 1:10){
    temp <- matchfun(string = "extra_se", obj = rawctl[,a])
    if(!is.na(temp)) { break() } }
 temp2<- as.numeric(rawctl[(temp+1):(temp+2),4])
 include.extra.se <- ifelse(sum(temp2)==0, FALSE, TRUE)
 # How many surveys have added se?
 n.extra.se = sum(temp2)

# Read Starter File and Change values
 starter.file <- SS_readstarter("starter.ss", verbose = FALSE)
 starter.file$detailed_age_structure <- 0
 SS_writestarter(starter.file,file="starter.ss",overwrite=T, verbose=FALSE, warn=FALSE)


 # Set up storage matrix
 Quant.out  <-define_matrix(N = Niter, ofl.yrs, depl.yr, n.survey, n.extra.se) 
 
 print("Getting Intial Sample")
 start.time <- Sys.time()
 # Do the initial SS runs
 for (i in 1:Niter)
 {
    changeM(ctl = control.name, para=parm.vec[i,1:2])
    changeH(ctl = control.name, para=parm.vec[i,3])
    changeDepl(dat = dat.name, para=parm.vec[i,4])  

    # Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1")}      
    
    rep.new       <- readLines(file.name)
    Quant.out[i,] <- getQuant(rep.new, parm=parm.vec[i,], ofl.yrs, depl.yr, 
                              n.extra.se, n.survey)
    
    #Rename the report file by rep number and move to a file to save for future needs
    move.files.fxn(rep.folder = rep.folder, sim.num=i)
      
    quant.list[[1]]   <- Quant.out
    save(quant.list, file=quant.file)
 }
 end.time <- Sys.time()
 print("Intial Sample Completed Taking Total Time of:")
 print(end.time - start.time)
 
 
 #Remove runs where depletion was not met
 find <- Quant.out[,"MissDep"] ==0 & Quant.out[,"Crash"] == 0
 Quant.out.good = Quant.out[find,]
  
 init.parms.list[[1]]   <- parm.vec #Save the initial parameter before removing them
 parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
 colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
 
 #Get the likelihood value for each run
 #perhaps change this to exp(-LLsurvey) instead or exp(-NLL) for total likelihood value
 likelihood <- exp(-Quant.out.good[,"NLL"])
 
 #Calculate the prior probability based on the initial good parameter vectors
 prior.all     <- get.prior(new.para=parm.vec, m.in, h.in, depl.in) 

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
          posterior.all    <- get.prior(new.para=parm.vec, m.in, h.in, depl.in)
          p                <- likelihood/sum(likelihood) 
          sample.wghts     <- p 
          sir.out          <- do.sir(Ncount=floor(0.25*Niter), input= parm.vec, wghts=sample.wghts)
          #Check for how many unique SIR draws 
          unq.draw         <- length(unique(sir.out$get.samp))         
      }
    
    #Here after sample weights are = likelihood*prior/pr
    if (ais > 1 )
      {
         ais.parm.vec       <- parm.vec
         p                  <- likelihood/sum(likelihood)
         posterior.all      <- get.new.posterior(new.para=ais.parm.vec, degree=5, m.in, h.in, depl.in)
         posterior          <- posterior.all$posterior
         prior.all          <- get.prior(new.para=ais.parm.vec, m.in, h.in, depl.in) 
         prior              <- prior.all$prior
         sample.wghts       <- get.new.wghts(p, prior, posterior)
         sir.out            <- do.sir(Ncount=floor(0.25*Niter),input= ais.parm.vec, wghts=sample.wghts)       
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
        
    print(c("Entropy:", round(entropy,4)))
    print(c("Effective N:", round(effN)))
    print(c("Expected Unique Points:", round(expN)))
    print(c("Max Weight:", round(maxW,4)))
    print(c("Variance of the rescaled:", varW))    
    if (entropy >= entropy.level) break()
                         
    #This is where sample from the new parameter values
    new.dists           <- fit.mvt(Niter, para=sir.vec, degree=5, m.in, depl.in)
    ais.parm.vec        <- do.call("cbind", new.dists)
    
    Quant.out  <-define_matrix(N = Niter, ofl.yrs, depl.yr, n.survey, n.extra.se) 
    #Do the initial SS runs
    for (i in 1:Niter)
    {
        changeM(ctl = control.name, para=ais.parm.vec[i, c("M.f", "M.m")])
        changeH(ctl = control.name, para=ais.parm.vec[i, "h"])
        changeDepl(dat = dat.name, para=ais.parm.vec[i, "depl"])
        
        #Run Simple Stock Synthesis
        if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
        if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1")}

        rep.new             <- readLines("ss_summary.sso")
        Quant.out[i,]       <- getQuant(rep.new, parm = ais.parm.vec[i,], ofl.yrs, depl.yr, n.extra.se, n.survey)
        quant.list[[ais+1]] <- Quant.out
        save(quant.list, file = quant.file)  
    }
    
    #Remove runs Crashed Runs and where Depletion was not Met
    find <- Quant.out[,"MissDep"] ==0 & Quant.out[,"Crash"] == 0
    Quant.out.good = Quant.out[find,]
 
    parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
    colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
    
    #Get the likelihood value for each run
    #trans.like <- exp(-(Quant.out.good$LL_survey-min(Quant.out.good$LL_survey)))
    likelihood <- exp(-Quant.out.good[,"NLL"])

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
 final.parm.vec      <- do.call("cbind", final.sir.vec)

 
 #Create Storage matrixes
 SB <- TotBio <- SmryBio <- SPR <- Exploitation <- as.data.frame(matrix(NA,nrow=length(all.yrs),ncol=final.Niter))
 colnames(SB) = colnames(TotBio) = colnames(SmryBio) = colnames(SPR) = colnames(Exploitation) = 1:final.Niter
 rownames(SB) = rownames(TotBio) = rownames(SmryBio) = rownames(SPR) = rownames(Exploitation) = all.yrs

 Bratio <-as.data.frame(matrix(NA,nrow=(length(all.yrs)-1),ncol=final.Niter))
 colnames(Bratio) = 1:final.Niter
 rownames(Bratio) = all.yrs[2]:all.yrs[length(all.yrs)]

 OFL <- ABC <- as.data.frame(matrix(NA,nrow=length(ofl.yrs),ncol=final.Niter))
 colnames(OFL) = colnames(ABC) = 1:final.Niter 
 rownames(OFL) = rownames(ABC) = ofl.yrs

 Survey <- matrix(NA, nrow = length(survey.yrs), ncol = final.Niter)
 
 #Do the final SS run
 Quant.out  <-define_matrix(N = final.Niter, ofl.yrs, depl.yr, n.survey, n.extra.se) 
 for (i in 1:final.Niter)
 {
    changeM(ctl = control.name, para=final.parm.vec[i,c("M.f", "M.m")])
    changeH(ctl = control.name, para=final.parm.vec[i,"h"])
    changeDepl(dat = dat.name, para=final.parm.vec[i,"depl"])
        
    #Run Simple Stock Synthesis
    if (tantalus == TRUE) { system("./SS -nohess > out.txt 2>&1")  }
    if (tantalus == FALSE){ shell("ss.exe -nohess > out.txt 2>&1") }
    
    rep.new       <- readLines(file.name)
    Quant.out[i,] <- getQuant(rep.new, parm=final.parm.vec[i,], ofl.yrs, depl.yr, n.extra.se, n.survey)
    quant.list[[Counter+2]] <- Quant.out
    save(quant.list, file=quant.file) 
       
    Summary   <- RepSumFxn(rep.new, n=i, all.yrs, hist.yrs, ofl.yrs)
    SB[,i]       <- Summary$SB
    SmryBio[,i]  <- Summary$SmryBio
    SPR[,i]      <- Summary$SPR
    TotBio[,i]   <- Summary$TotBio
    Bratio[,i]   <- Summary$Bratio
    OFL[,i]      <- Summary$OFL
    ABC[,i]      <- Summary$ForeCat
    Exploitation[,i] <- c(Catch, Summary$ForeCat) / Summary$TotBio
    Survey[,i]   <- as.numeric(Summary$Survey)
     
    #Rename the report file by rep number and move to a file to save for future needs
    move.files.fxn(rep.folder = rep.folder, sim.num=i) 
    
 }
 end.time <- Sys.time()
 print(end.time - start.time)
 
 rownames(Survey) = names(Summary$Survey)
 #Remove runs Crashed Runs and where Depletion was not met
 index <- Quant.out[,"MissDep"] ==0 & Quant.out[,"Crash"] == 0
 
 Quant.out.good <- Quant.out[index,]
 rep.list[[1]]  <- TotBio[,index]
 rep.list[[2]]  <- SB[,index]
 rep.list[[3]]  <- SmryBio[,index]
 rep.list[[4]]  <- Bratio[,index]
 rep.list[[5]]  <- SPR[,index]
 rep.list[[6]]  <- OFL[,index]
 rep.list[[7]]  <- ABC[,index]
 rep.list[[8]]  <- Exploitation[,index]
 rep.list[[9]]  <- Catch
 rep.list[[10]] <- Survey
 
 names(rep.list) <- c("TotBio", "SB", "SmryBio", "Bratio", "SPR", "OFL", "ABC", "Exploitation", "Catch", "Survey")
 
 parm.vec <- as.data.frame(Quant.out.good[,1:4]) #Replace the parm vector with only the good runs
 colnames(parm.vec) <- c("M.f", "M.m", "h", "depl")
 
 #Get the likelihood value for each run
 likelihood    <- exp(-Quant.out.good[,"NLL"])
 prior.all     <- get.prior(new.para=parm.vec[,1:4], m.in, h.in, depl.in)  
 posterior.all <- get.new.posterior(new.para=parm.vec[,1:4], degree=5, m.in, h.in, depl.in)
 
 #Save the final trajectories, liklihood values, and probabilities
 posterior.prior.list[[Counter+2]]  <- c(prior.all, posterior.all)
 quant.good.list[[Counter+2]]   <- Quant.out.good
 parm.list[[Counter+2]] <- parm.vec

 #Write saved files with year index
 save(quant.good.list, file=quant.good.file)
 save(rep.list, file= rep.file)
 save(posterior.prior.list, file=posterior.prior.file)
 save(parm.list, file=parameters)
 
 pdf(paste0(save.folder,"/posteriors.pdf"),width=7,height=7,)
 par(mfrow=c(2,2),oma=c(3,3,4,3))
 hist(parm.vec[,1],xlab= paste("Natural Mortality F (mean=", round(mean(parm.vec[,1]),2),"sd=",round(sd(parm.vec[,1]),2),")", sep=" "), main="")
 hist(parm.vec[,2],xlab= paste("Natural Mortality M (mean=",round(mean(parm.vec[,2]),2),"sd=",round(sd(parm.vec[,2]),2),")", sep=" "),main="")
 hist(parm.vec[,3],xlab= paste("Steepness (mean=",round(mean(parm.vec[,3]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 hist(parm.vec[,4],xlab= paste("Depletion Target (,",depl.yr," (mean=",round(mean(parm.vec[,4]),2),"sd=",round(sd(parm.vec[,4]),2),")",sep=" "),main="") 
 mtext("Posterior Distributions", side = 3, outer=T)
 dev.off()

 create.Plots(dir = save.folder, rep.list, parm.list, quant.list = quant.good.list, dat.name,
                all.yrs, ofl.yrs, hist.yrs, depl.in, m.in, h.in, n.extra.se, n.survey)
 
 return("XSSS Completed")
}
