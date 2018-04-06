##' Changes the natural mortality value for female and male fish in the control file 
##' @param ctl
##' @param para M vector for female and male values for the nth simulation
##' @author Chantel Wetzel
##' @export
 changeM <- function(ctl, para) 
 {
    para <- round(para,4)
    #ctl.new <- readLines(ctl)
    #M.line  <- strsplit(ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)],"[[:blank:]]+")[[1]]
    #M.line[3] <- para[1]
    #ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    #
    #M.line  <- strsplit(ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)],"[[:blank:]]+")[[1]]
    #M.line[3] <- para[2]
    #ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    #write(ctl.new, ctl)

    SS_changepars(ctlfile=ctl, newctlfile=ctl, strings=c("NatM_p_1_Fem_GP_1","NatM_p_1_Mal_GP_1") , 
        newvals=c(para[1], para[2]), estimate=FALSE, verbose=FALSE)
 }