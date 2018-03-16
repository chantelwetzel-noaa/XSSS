##' Seed Function - this function will generate a seed if one is not read from a file
##' @author Chantel Wetzel
##' @export 
  get.seed <- function(){ 
    Sys.sleep(1)
    t <- as.numeric(Sys.time()) ; seed <- (t-floor(t))*1e8 # set seed base upon epochs
    return(seed)
  }