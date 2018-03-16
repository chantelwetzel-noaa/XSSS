##' Function that moves and saves the needed files
##' sim.num simulation number that is appended to the report file
##' @author Chantel Wetzel
##' @export
 move.files.fxn <- function(sim.num)
 {
    #Rename the report file by rep number and move to a file to save for future needs
    file.rename(paste(directory,"Report.sso",sep=""),paste(directory,"Report",sim.num,".sso",sep=""))
    file.copy(paste(directory,"Report",sim.num,".sso",sep=""), paste(rep.folder,"/Report",sim.num,".sso",sep=""),overwrite=T)
    file.remove(paste(directory,"Report",sim.num,".sso",sep=""))
 }