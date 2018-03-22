##' Function that moves and saves the needed files
##' sim.num simulation number that is appended to the report file
##' @author Chantel Wetzel
##' @export
 move.files.fxn <- function(rep.folder, sim.num)
 {
 	#name = ifelse(SS_versionNumeric >= 3.3, "ss_summary.sso", "Report.sso")

    #Rename the report file by rep number and move to a file to save for future needs
    #file.rename(paste(directory,"ss_summary.sso",sep=""),paste(directory,"ss_summary",sim.num,".sso",sep=""))
    #file.copy(paste(directory,  "ss_summary",sim.num,".sso",sep=""), paste(rep.folder,"/ss_summary",sim.num,".sso",sep=""),overwrite=T)
    #file.remove(paste(directory,"ss_summary",sim.num,".sso",sep=""))

    file.rename("ss_summary.sso",paste0("ss_summary",sim.num,".sso"))
    file.copy  (paste0("ss_summary",sim.num,".sso"), paste0(rep.folder,"/ss_summary",sim.num,".sso"),overwrite=T)
    file.remove(paste0("ss_summary",sim.num,".sso"))
 }