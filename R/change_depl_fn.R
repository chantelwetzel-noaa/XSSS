##' Changes the depletion value for the depletion year in the data file 
##' @param dat 
##' @param para depletion parameter value for the nth simulation
##' @author Chantel Wetzel
##' @export
 changeDepl <- function(dat, para) 
 {
    para   <- round(para,4)
    dat.new    <- readLines(dat)
    depl.line <- strsplit(dat.new[grep("FinalDepl",dat.new)]," ")[[1]]
    depl.line[4] <- para
    dat.new[grep("FinalDepl",dat.new)][[1]] <- paste(depl.line,collapse=" ")
    write(dat.new, dat)
    
 }