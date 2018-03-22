##' 
##' @param N
##' @param ofl.yrs
##' @param depl.yr
##' @author Chantel Wetzel
##' @export 
define_matrix <- function(N, ofl.yrs, depl.yr){
	names <- c("M_f","M_m","h","depl","LnR0","extra_se","SB0",paste0("SPB_",ofl.yrs[1]),paste0("Dep_TgtYr_",depl.yr),
              paste0("Dep_", ofl.yrs[1]), paste0("OFL_",ofl.yrs[2]), paste0("AdjCatch_",ofl.yrs[2]),
              "SBMSY/SB0","SSB_Unfished","SmryBio_Unfished","SSB_SPRtgt","Fstd_SPRtgt",
              "TotYield_SPRtgt","SSB_MSY","Fstd_MSY","NLL","LL_survey","Crash","MissDep", "Survey_Q")
 	out  <-as.matrix(matrix(NA, nrow=N, ncol=length(names)))
 	colnames(out) <- names
 	return(out)
}