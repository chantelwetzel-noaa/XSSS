##' Function to pull quantities
##' @param string
##' @param obj
##' @param substr1
##' @author Ian Taylor
##' @export
matchfun <- function(string, obj, substr1=TRUE)
{
	# return a line number from the report file (or other file)
	# sstr controls whether to compare subsets or the whole line
	match(string, if(substr1){substring(obj,1,nchar(string))}else{obj} )
}