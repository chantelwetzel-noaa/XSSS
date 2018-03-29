##' Get Final Run Model Output   
##' @param n simulation number
##' @param rep.new report file name to pull summary values from
##' @param all.yrs 
##' @param hist.yrs 
##' @param ofl.yrs 
##' @author Chantel Wetzel
##' @export
RepSumFxn<- function(rep.new, n, all.yrs, hist.yrs, ofl.yrs)
 {
     
    #if (ssver < 3.30){
    #    TotBio= mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste(1, x,"TIME",sep=" "),rep.new)]," ")[[1]][5]),x=hist.yr)
    #    SB = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste("SPB_",x,sep=""),rep.new)]," ")[[1]][3]),x=all.yrs) 
    #    Bratio= mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste("Bratio_",x,sep=""),rep.new)], " ")[[1]][3]), x=all.yrs[2]:all.yrs[length(all.yrs)])
    #    OFL = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][3]), x=ofl.yrs)
    #    ForeCat= mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][3]), x=ofl.yrs)    
    #}

    #if (ssver >= 3.30){
    TotBio  = mapply(function(x) TotBio = as.numeric(strsplit(rep.new[grep(paste0("TotBio_",x),rep.new)]," ")[[1]][2]),x=all.yrs)
    SB      = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste0("SSB_",x),rep.new)]," ")[[1]][2]),x=all.yrs) 
    SmryBio = mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste0("SmryBio_",x),rep.new)]," ")[[1]][2]),x=all.yrs) 
    SPR     = 1 - mapply(function(x) SB = as.numeric(strsplit(rep.new[grep(paste0("SPRratio_",x),rep.new)]," ")[[1]][2]),x=all.yrs) 
    Bratio  = mapply(function(x) Bratio = as.numeric(strsplit(rep.new[grep(paste0("Bratio_",x),rep.new)], " ")[[1]][2]), x=all.yrs[2]:all.yrs[length(all.yrs)])
    OFL     = mapply(function(x) OFL = as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",x,sep=""),rep.new)], " ")[[1]][2]), x=ofl.yrs)
    ForeCat = mapply(function(x) ForeCat = as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",x,sep=""),rep.new)], " ")[[1]][2]), x=ofl.yrs)    
    #}

    start = grep("#_survey_stdev", rep.new) + 1
    end   = grep("#_Biomass", rep.new) - 1

    df <- rep.new[start:end]
    df <- strsplit(df, "[[:blank:]]+") ## Split by whitespace and collapse (+)
    df <- as.list(df) 
    df <- do.call("rbind", df)
    Survey = noquote(df[,5]) 
    names(Survey) = noquote(df[,1])   

    Summary <- list()
    Summary[[1]] <- TotBio
    Summary[[2]] <- SB
    Summary[[3]] <- Bratio
    Summary[[4]] <- OFL
    Summary[[5]] <- ForeCat
    Summary[[6]] <- SPR
    Summary[[7]] <- SmryBio
    Summary[[8]] <- Survey

    names(Summary) <- c("TotBio", "SB","Bratio","OFL", "ForeCatch", "SPR", "SmryBio", "Survey")
    return(Summary)
 }