##' Create the quant matrix
##' @param N
##' @param ofl.yrs
##' @param depl.yr
##' @param n.extra.se
##' @author Chantel Wetzel
##' @export 
define_matrix <- function(N, ofl.yrs, depl.yr, n.extra.se){

	names <- c("M_f","M_m","h","depl","LnR0","extra_se","SB0",paste0("SPB_",ofl.yrs[1]),paste0("Dep_TgtYr_",depl.yr),
              paste0("Dep_", ofl.yrs[1]), paste0("OFL_",ofl.yrs[2]), paste0("AdjCatch_",ofl.yrs[2]),
              "SBMSY/SB0","SSB_Unfished","SmryBio_Unfished","SSB_SPRtgt","Fstd_SPRtgt",
              "TotYield_SPRtgt","SSB_MSY","Fstd_MSY","NLL","LL_survey","Crash","MissDep", 
              rep(paste0("Survey_Q_",1:n.extra.se), n.extra.se))
 	out  <-as.matrix(matrix(NA, nrow=N, ncol=length(names)))
 	colnames(out) <- names
 	return(out)
}