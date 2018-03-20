
##' Function List
##; @author Chantel Wetzel
##' @export
dir = paste0(file.path, "/R/")

source(paste0(dir, "rbeta_ab.R"))
source(paste0(dir, "pars_truncbeta.R"))
source(paste0(dir, "do_sir.R"))
source(paste0(dir, "get_prior.R"))
source(paste0(dir, "get_new_posterior.R"))
source(paste0(dir, "get_new_wghts.R"))
source(paste0(dir, "fit_mvt.R"))
source(paste0(dir, "change_M.R"))
source(paste0(dir, "change_H.R"))
source(paste0(dir, "change_Depl.R"))
source(paste0(dir, "get_Quant.R"))
source(paste0(dir, "Rep_Summary_Fxn.R"))
source(paste0(dir, "move_file_fxn.R"))
source(paste0(dir, "get_seed.R"))
source(paste0(dir, "Match_Fn.R"))