##' Function to pull quantities
##' @param rep.new
##' @param parm
##' @param ofl.yrs
##' @param depl.yrs
##' @param n.extra.se
##' @param n.survey
##' @author Chantel Wetzel
##' @export
 getQuant <- function(rep.new, parm, ofl.yrs, depl.yr, n.extra.se, n.survey)
 {
     
                                        
    sb.name  <- "SSB"
    unfished <- "unfished"
    index    <- 2
    spr.tgt  <- "SPR"
    yield    <- "Dead_Catch"

    se = q = numeric(0)
    for (i in 1:n.survey){
        q[i] = as.numeric(strsplit(rep.new[grep("LnQ_base",rep.new)], " ")[[i]][index])
    }        

    if (n.extra.se > 0){
        for(i in 1:n.extra.se){
            se[i] = as.numeric(strsplit(rep.new[grep("Q_extraSD",rep.new)], " ")[[i]][index])
        }
    }
    if(n.extra.se == 0 ) {se = NA}

    out <- list() 
    out[1]  <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][index])
    out[2]  <- as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][index])
    out[3]  <- as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][index])
    out[4]  <- parm["depl"]
    out[5]  <- as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][index])
    out[6]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[7]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])
    out[8]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",depl.yr),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[9]  <- as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_",ofl.yrs[1]),rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name, "_Initial"),rep.new)], " ")[[1]][index])
    out[10] <- as.numeric(strsplit(rep.new[grep(paste0("OFLCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out[11] <- as.numeric(strsplit(rep.new[grep(paste0("ForeCatch_",ofl.yrs[2]),rep.new)], " ")[[1]][index])
    out[12] <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])/as.numeric(strsplit(rep.new[grep(paste0(sb.name,"_Initial"),rep.new)], " ")[[1]][index])
    out[13] <- as.numeric(strsplit(rep.new[grep(paste0("SSB_",unfished), rep.new)], " ")[[1]][index])
    out[14] <- as.numeric(strsplit(rep.new[grep(paste0("SmryBio_",unfished),rep.new)], " ")[[1]][index])
    out[15] <- as.numeric(strsplit(rep.new[grep(paste0("SSB_", spr.tgt),rep.new)], " ")[[1]][index])
    out[16] <- as.numeric(strsplit(rep.new[grep(paste0("Fstd_", spr.tgt),rep.new)], " ")[[1]][index])
    out[17] <- as.numeric(strsplit(rep.new[grep(paste0(yield, "_", spr.tgt),rep.new)], " ")[[1]][index])
    out[18] <- as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][index])
    out[19] <- as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][index])
    out[20] <- as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    out[21] <- as.numeric(strsplit(rep.new[grep("Survey",rep.new)], " ")[[1]][2])
    out[22] <- as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
    out[23] <- abs(as.numeric(out[4]) - as.numeric(out[8])) > 0.01
    ind = length(out) + 1
    for(i in ind:(ind+n.survey-1))  { out[i] <- q[i-ind + 1] }
    ind = length(out) + 1
    for(i in ind:(ind+n.extra.se-1)){ out[i] <- se[i - ind + 1] } 
    
    
    out = do.call("cbind", out)
    out = unlist(out)
    return(out)
 }