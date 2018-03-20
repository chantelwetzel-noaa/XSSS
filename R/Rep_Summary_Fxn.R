##' Get Final Run Model Output   
##' @param n simulation number
##' @param rep.new report file name to pull summary values from
##' @author Chantel Wetzel
##' @export
RepSumFxn<- function(rep.new, n)
 {
    tot.yrs <- c(hist.yrs, ofl.yrs)
          
    TotBio= mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste(1, x,"TIME",sep=" "),rep.new)]," ")[[1]][5]),x=hist.yr)
    SB = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste("SPB_",x,sep=""),rep.new)]," ")[[1]][3]),x=tot.yrs) 
    Bratio= mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste("Bratio_",x,sep=""),rep.new)], " ")[[1]][3]),x=tot.yrs[2]:tot.yrs[length(tot.yrs)])
    OFL = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][3]),x=ofl.yrs)
    ForeCat= mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][3]),x=ofl.yrs)    

    RepSummary <- list()
    RepSummary[[1]] <- TotBio
    RepSummary[[2]] <- SB
    RepSummary[[3]] <- Bratio
    RepSummary[[4]] <- OFL
    RepSummary[[5]] <- ForeCat

    names(RepSummary) <- c("TotBio", "SB","Bratio","OFL", "ForeCatch")
    return(RepSummary)
 }