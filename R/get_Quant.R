##' This function is currently not being used but I have kept it in here just in case I want to revert back 
##' @param n simulation number
##' @param parm
##' @param N
##' @author Chantel Wetzel
##' @export
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