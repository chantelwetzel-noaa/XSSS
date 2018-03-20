##' Function to pull quantities
##' @param n simulation number
##' @param parm
##' @param N
##' @param ssver model version
##' @author Chantel Wetzel
##' @export
 getQuant <- function(rep.new,n,parm,N,ssver)
 {
    if (n == 1) {
    get     <- c("M_f","M_m","h","depl","LnR0","SB0",paste0("SPB_",ofl.yrs[1]),paste0("Dep_TgtYr_",depl.yr),
                  paste0("Dep_", ofl.yrs[1]), paste0("OFL_",ofl.yrs[2]), paste0("AdjCatch_",ofl.yrs[2]),
                  "SBMSY/SB0","SSB_Unfished","SmryBio_Unfished","SSB_SPRtgt","Fstd_SPRtgt",
                  "TotYield_SPRtgt","SSB_MSY","Fstd_MSY","NLL","LL_survey","Crash","R0_init","MissDep")
    out  <-as.data.frame(matrix(NA,nrow=N,ncol=length(get)))
    colnames(out) <- get
    }   
                                        
    sb.name  <- ifelse(ssvers >=3.3, "SSB", "SPB")
    unfished <- ifelse(ssvers >=3.3, "unfished", "Unfished")
    index    <- ifelse(ssvers >=3.3, 2, 3)
    spr.tgt  <- ifelse(ssvers >=3.3, "SPR", "SPRtgt")
    yield    <- ifelse(ssvers >=3.3, "Dead_Catch", "TotYield")


    out[n,1]  <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][index])
    out[n,2]  <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][index])
    out[n,3]  <- as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][index])
    out[n,4]  <- parm[,4]
    out[n,5]  <- as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][index])
    out[n,6]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[n,7]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])
    out[n,8]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",depl.yr),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[n,9]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name, "_Initial"),rep.new)], " ")[[1]][index])
    out[n,10] <- as.numeric(strsplit(rep.new[grep(paste0("OFLCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out[n,11] <- as.numeric(strsplit(rep.new[grep(paste0("ForeCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out[n,12] <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[n,13] <- as.numeric(strsplit(rep.new[grep(paste0("SSB_",unfished), rep.new)], " ")[[1]][index])
    out[n,14] <- as.numeric(strsplit(rep.new[grep(paste0("SmryBio_",unfished),rep.new)], " ")[[1]][index])
    out[n,15] <- as.numeric(strsplit(rep.new[grep(paste0("SSB_", spr.tgt),rep.new)], " ")[[1]][index])
    out[n,16] <- as.numeric(strsplit(rep.new[grep(paste0("Fstd_", spr.tgt),rep.new)], " ")[[1]][index])
    out[n,17] <- as.numeric(strsplit(rep.new[grep(paste0(yield, "_", spr.tgt),rep.new)], " ")[[1]][index])
    out[n,18] <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])
    out[n,19] <- as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][index])
    out[n,20] <- as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    out[n,21] <- as.numeric(strsplit(rep.new[grep("Survey",rep.new)], " ")[[1]][2])
    out[n,22] <- as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
    out[n,23] <- as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8]) #This needs to be initial R0
    out[n,24] <- abs(out[n,4] - out[n,8]) > 0.01
    
    return(out)
 }