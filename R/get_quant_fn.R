##' Function to pull quantities
##' @param parm
##' @param ssver model version
##' @author Chantel Wetzel
##' @export
 getQuant <- function(rep.new, parm, ofl.yrs, depl.yr, read.se, ssver)
 {
    out <- list()  
                                        
    sb.name  <- ifelse(ssver >=3.3, "SSB", "SPB")
    unfished <- ifelse(ssver >=3.3, "unfished", "Unfished")
    index    <- ifelse(ssver >=3.3, 2, 3)
    spr.tgt  <- ifelse(ssver >=3.3, "SPR", "SPRtgt")
    yield    <- ifelse(ssver >=3.3, "Dead_Catch", "TotYield")

    se = NA
    if(read.se){ se = as.numeric(strsplit(rep.new[grep("Q_extraSD",rep.new)], " ")[[1]][index]) }


    out$M.f <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][index])
    out$M.m <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][index])
    out$h <- as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][index])
    out$depl <- parm["depl"]
    out$LnR0 <- as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][index])
    out$SurveyExtraSE <- se
    out$SB0 <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out$SBfinal <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])
    out$Depltgtyr <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",depl.yr),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out$Deplfinal  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name, "_Initial"),rep.new)], " ")[[1]][index])
    out$ofl <- as.numeric(strsplit(rep.new[grep(paste0("OFLCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out$forecatch <- as.numeric(strsplit(rep.new[grep(paste0("ForeCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out$SBmsySB0 <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out$SBUnfished <- as.numeric(strsplit(rep.new[grep(paste0("SSB_",unfished), rep.new)], " ")[[1]][index])
    out$SumBioUnfished <- as.numeric(strsplit(rep.new[grep(paste0("SmryBio_",unfished),rep.new)], " ")[[1]][index])
    out$SBSPRtgt <- as.numeric(strsplit(rep.new[grep(paste0("SSB_", spr.tgt),rep.new)], " ")[[1]][index])
    out$FSPRtgt <- as.numeric(strsplit(rep.new[grep(paste0("Fstd_", spr.tgt),rep.new)], " ")[[1]][index])
    out$YeildSPRtgt <- as.numeric(strsplit(rep.new[grep(paste0(yield, "_", spr.tgt),rep.new)], " ")[[1]][index])
    out$SBmsy <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])
    out$Fmsy <- as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][index])
    out$Like <- as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    out$SurveyLike <- as.numeric(strsplit(rep.new[grep("Survey",rep.new)], " ")[[1]][2])
    out$CrashPen <- as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
    out$HitDepl <- abs(as.numeric(out$depl) - as.numeric(out$Depltgtyr)) > 0.01
    out$LnQ <- as.numeric(strsplit(rep.new[grep("LnQ_base",rep.new)], " ")[[1]][index])
    
    out = do.call("cbind", out)
    out = unlist(out)
    return(out)
 }