##' Get Final Run Model Output   
##' @param n simulation number
##' @param rep.new report file name to pull summary values from
##' @param all.yrs 
##' @param hist.yrs 
##' @param ofl.yrs 
##' @param ssvers 
##' @author Chantel Wetzel
##' @export
RepSumFxn<- function(rep.new, n, all.yrs, hist.yrs, ofl.yrs, ssver)
 {
    #tot.yrs <- c(hist.yrs, ofl.yrs)
     
    if (ssver < 3.30){
        TotBio= mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste(1, x,"TIME",sep=" "),rep.new)]," ")[[1]][5]),x=hist.yr)
        SB = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste("SPB_",x,sep=""),rep.new)]," ")[[1]][3]),x=all.yrs) 
        Bratio= mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste("Bratio_",x,sep=""),rep.new)], " ")[[1]][3]), x=all.yrs[2]:all.yrs[length(all.yrs)])
        OFL = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][3]), x=ofl.yrs)
        ForeCat= mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][3]), x=ofl.yrs)    
    }

    if (ssver >= 3.30){
        TotBio  = mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste0("TotBio_",x),rep.new)]," ")[[1]][2]),x=all.yrs)
        SB      = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste0("SSB_",x),rep.new)]," ")[[1]][2]),x=all.yrs) 
        Bratio  = mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste0("Bratio_",x),rep.new)], " ")[[1]][2]), x=all.yrs[2]:all.yrs[length(all.yrs)])
        OFL     = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][2]), x=ofl.yrs)
        ForeCat = mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][2]), x=ofl.yrs)    
    }       

    RepSummary <- list()
    RepSummary[[1]] <- TotBio
    RepSummary[[2]] <- SB
    RepSummary[[3]] <- Bratio
    RepSummary[[4]] <- OFL
    RepSummary[[5]] <- ForeCat

    names(RepSummary) <- c("TotBio", "SB","Bratio","OFL", "ForeCatch")
    return(RepSummary)
 }