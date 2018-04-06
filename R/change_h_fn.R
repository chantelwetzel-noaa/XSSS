##' Changes the steepness value in the control file 
##' @param para steepness parameter value for the nth simulation
##' @author Chantel Wetzel
##' @export 
changeH <- function(ctl, para)
 {
    para <- round(para,4)
    #ctl.new <- readLines(ctl)
    #h.line  <- strsplit(ctl.new[grep("_steep",ctl.new)],"[[:blank:]]+")[[1]]
    #h.line[3] <- para
    #ctl.new[grep("_steep",ctl.new)] <- paste(h.line,collapse=" ")
    #write(ctl.new, ctl) 

    SS_changepars(ctlfile=ctl, newctlfile=ctl, strings=c("_steepn") , 
        newvals=para, estimate=FALSE, verbose=FALSE)  
 }