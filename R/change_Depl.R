##' Changes the depletion value for the depletion year in the data file 
##' @param para depletion parameter value for the nth simulation
##' @author Chantel Wetzel
##' @export
 changeDepl <- function(para) 
 {
    para   <- round(para,4)
    dat    <- readLines(paste0(directory, dat.name))
    depl.line <- strsplit(dat[grep("FinalDepl",dat)]," ")[[1]]
    depl.line[4] <- para
    dat[grep("FinalDepl",dat)][[1]] <- paste(depl.line,collapse=" ")
    write(dat,paste(directory,dat.name,sep=""))
    
 }