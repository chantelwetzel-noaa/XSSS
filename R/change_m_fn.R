##' Changes the natural mortality value for female and male fish in the control file 
##' @param para M vector for female and male values for the nth simulation
##' @author Chantel Wetzel
##' @export
 changeM <- function(ctl, para) 
 {
    para <- round(para,4)
    ctl.new <- readLines(ctl)
    M.line  <- strsplit(ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)],"[[:blank:]]+")[[1]]
    M.line[c(3,4)] <- para[1]
    ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    M.line  <- strsplit(ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)],"[[:blank:]]+")[[1]]
    M.line[c(3,4)] <- para[2]
    ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)] <- paste(M.line,collapse=" ")
    write(ctl.new, ctl)
 }