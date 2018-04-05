##' Changes the steepness value in the control file 
##' @param para steepness parameter value for the nth simulation
##' @author Chantel Wetzel
##' @export 
changeH <- function(ctl, para)
 {
    para <- round(para,4)
    ctl.new <- readLines(ctl)
    h.line  <- strsplit(ctl.new[grep("BH_steep",ctl.new)],"[[:blank:]]+")[[1]]
    h.line[c(3,4)] <- para
    ctl.new[grep("_steep",ctl.new)] <- paste(h.line,collapse=" ")
    write(ctl.new, ctl)   
 }