##' 
##' @param Ncount
##' @param input
##' @param wghts
##' @author Chantel Wetzel
##' @export 
do.sir <- function(Ncount, input, wghts)
 {
    
    get.samp   <- sample(1:length(wghts),size=Ncount,replace=T,prob=wghts)
    sir.vec    <- input[get.samp,]
      
    output            <- NULL
    output$sir.vec    <- sir.vec
    output$get.samp   <- get.samp
    return(output)
 } 