##' Changes the steepness value in the control file 
##' @param para steepness parameter value for the nth simulation
##' @author Chantel Wetzel
##' @export 
changeH <- function(para)
 {
    para <- round(para,4)
    ctl.new <- readLines(paste(directory,control.name,sep=""))
    h.line  <- strsplit(ctl.new[grep("SR_steep",ctl.new)]," ")[[1]]
    h.line[c(3,4)] <- para
    ctl.new[grep("SR_steep",ctl.new)] <- paste(h.line,collapse=" ")
    write(ctl.new,paste(directory,control.name,sep=""))   
 }