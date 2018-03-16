##' Calculate the new sampling weights based on the likelihood, prior, and posterior values
##' @param like
##' @param prior
##' @param posterior
##' @author Chantel Wetzel
##' @export 
  get.new.wghts <- function(like,prior,posterior)
  {
    sample.wghts <- (like*prior)/posterior
    new.p        <- sample.wghts/sum(sample.wghts)
    return(new.p)
  }