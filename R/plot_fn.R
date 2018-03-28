##' Changes the natural mortality value for female and male fish in the control file 
##' @param dir the folder where plots and tables will be saved
##' @param rep.list timeseries list
##' @param parm.list
##' @param quant.list
##' @param all.yrs
##' @param ofl.yrs
##' @param hist.yrs
##' @param depl.in
##' @param m.in
##' @param h.in
##' @param n.extra.se
##' @param n.survey
##' @author Chantel Wetzel
##' @export
create.Plots <- function(dir = save.folder, rep.list, parm.list, quant.list,
						 all.yrs, ofl.yrs, hist.yrs, depl.in, m.in, h.in, n.extra.se, n.survey){

	#Get the colors for the polygons
	rc <- function(n,alpha) {
	    x <- seq(0,1,length=n)
	    r <- 1/(1+exp(20-35*x))
	    g <- pmin(pmax(0,-0.80+6*x-5*x^2),1)
	    b <- dnorm(x,0.25,0.15)/max(dnorm(x,0.25,0.15))
	    rgb.m <- matrix(c(r,g,b),ncol=3)
	    rich.vector <- apply(rgb.m,1,function(v) rgb(v[1],v[2],v[3],alpha=alpha))
	}
	shadecol <- rgb(0,0,0,alpha=0.10)

	# Create specialized plots
	pngfun <- function(file,w=7,h=7,pt=12){
	  file <- file.path(dir, file)
	  png(filename=file,
	      width=w,height=h,
	      units='in',res=300,pointsize=pt)
	}

	comma 		   <- function(x, digits=0) { formatC(x, big.mark=",", digits, format = "f") }

	# Compare prior and posterior distributions
	pngfun(file = "/Parameters.png")
	col.vec = c("red", "green", "blue")
	par(mfrow=c(2,2),mar=c(4,4,4,4), oma=c(1,1,1,1))
	val = quant.list; xx = length(val) 
	a <- density(parm.list[[1]]$M.f, bw=0.03) # Prior
	b <- density(val[[1]][,"M_f"],   bw=0.03) # Post Model Pre-Data
	c <- density(val[[xx]][,"M_f"],  bw=0.03) # Posterior
	xlim = max(a$x, b$x, c$x)
	plot(a, xlab="Female M", main='', ylim=c(0,max(a$y,b$y,c$y)), lwd=2, xlim=c(0,xlim), col = col.vec[1])
	lines(b, col= col.vec[2],lwd=2, lty = 3)
	lines(c,col=col.vec[3], lwd=2, lty=1)
	legend('topright',legend=c("Prior","Post-Model","Posterior"), col=col.vec,lwd=c(2,2,2),bty='n',lty=c(1,3,1))
	legend("right", legend = paste("Posterior Median = ", round(median(val[[xx]][,"M_f"]),3)) , bty = 'n')

	a <- density(parm.list[[1]]$M.m,bw=0.03) 
	b <- density(val[[1]][,"M_m"],  bw=0.03) 
	c <- density(val[[xx]][,"M_m"], bw=0.03) 
	xlim = max(a$x, b$x, c$x)
	plot(a, xlab="Male M", main='', ylim=c(0,max(a$y,b$y,c$y)), lwd=2, xlim=c(0,xlim), col = col.vec[1])
	lines(b, col= col.vec[2],lwd=2, lty = 3)
	lines(c,col=col.vec[3], lwd=2, lty=1)
	legend("right", legend = paste("Posterior Median = ", round(median(val[[xx]][,"M_m"]),3)) , bty = 'n')

	bw.val = 0.05
	a <- density(parm.list[[1]]$h,bw=bw.val, from=0.2, to=1)
	b <- density(val[[1]][,"h"],  bw=bw.val, from=0.2, to=1)
	c <- density(val[[xx]][,"h"],  bw=bw.val, from=0.2, to=1)
	xlim = max(a$x, b$x, c$x)
	plot(a,xlab="Steepness",main='',ylim=c(0, max(a$y,b$y,c$y)), lwd=2, xlim=c(0.2,1), col= col.vec[1])
	lines(b, col= col.vec[2],lwd=2, lty = 3)
	lines(c,col=col.vec[3], lwd=2, lty=1)
	legend("right", legend = paste("Posterior Median = ", round(median(val[[xx]][,"h"]),3)) , bty = 'n')
	
	a <- density(parm.list[[1]]$depl,from=0,to=1, bw=0.10)
	b <- density(val[[1]][,"depl"],  from=0,to=1, bw=0.10)
	c <- density(val[[xx]][,"depl"], from=0,to=1, bw=0.10)
	plot(a,xlab="Depletion in 2000",main='',ylim=c(0,max(a$y,b$y,c$y)),lwd=2,xlim=c(0,1), col= col.vec[1])
	lines(b, col= col.vec[2],lwd=2, lty = 3); lines(c,col=col.vec[3], lwd=2, lty=1)
	legend("right", legend = paste("Posterior Median = ", round(median(val[[xx]][,"depl"]),3)) , bty = 'n')
	dev.off()

	# Plot Parameter Correlation
	pngfun(file = "/Parameter_Correlation.png")
	par(mfrow=c(1,1),mar=c(4,4,4,4), oma=c(1,1,1,1))
	pairs(val[[xx]][,1:4])
	dev.off()

	# Total Spawning Output
	shadecol <- rc(n=2,alpha=0.10)
	linecol  <- rc(n=2,alpha=1)
	sb = apply(rep.list$SB, 1, quantile, c(0.025, 0.50, 0.975))
	pngfun (file = "/Spawning_Output.png")
	par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(2,2,2,2))
	plot(all.yrs, sb[2,], ylim=c(0, max(sb)), ylab="SB",xlab="Year",type='n',xaxs='i',axes=F)
	xx <- c(all.yrs, rev(all.yrs))
	yy <- c(sb[1,],rev(sb[3,]))
	polygon(xx, yy, col=shadecol[1], border='NA')
	lines(all.yrs, sb[2,],lty=1,col=linecol[1],lwd=2)
	box()
	axis(side=1); axis(side=2)
	legend("topright", legend = "Median with 95% Interval", bty = 'n')
	legend("bottomleft", legend = paste("Median(SB) ", hist.yrs[length(hist.yrs)], " = ", round(sb[2,length(hist.yrs)])) , bty = 'n')
	dev.off() 

	# Total Biomass
	shadecol <- rc(n=2,alpha=0.10)
	linecol  <- rc(n=2,alpha=1)
	sb = apply(rep.list$TotBio, 1, quantile, c(0.025, 0.50, 0.975))
	pngfun (file = "/Total_Biomass.png")
	par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(2,2,2,2))
	plot(all.yrs, sb[2,], ylim=c(0, max(sb)), ylab="Total Biomass",xlab="Year",type='n',xaxs='i',axes=F)
	xx <- c(all.yrs, rev(all.yrs))
	yy <- c(sb[1,],rev(sb[3,]))
	polygon(xx, yy, col=shadecol[1], border='NA')
	lines(all.yrs, sb[2,],lty=1,col=linecol[1],lwd=2)
	box()
	axis(side=1); axis(side=2)
	legend("topright", legend = "Median with 95% Interval", bty = 'n')
	legend("bottomleft", legend = paste("Median(Total Biomass) ", hist.yrs[length(hist.yrs)], " = ", round(sb[2,length(hist.yrs)])) , bty = 'n')
	dev.off() 

	# Depletion
	shadecol <- rc(n=2,alpha=0.10)
	linecol  <- rc(n=2,alpha=1)
	sb = apply(rep.list$Bratio, 1, quantile, c(0.025, 0.50, 0.975))
	pngfun (file = "/Depletion.png")
	par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(2,2,2,2))
	yr = all.yrs[2:length(all.yrs)]
	plot(yr, sb[2,], ylim=c(0, 1), ylab="Depletion",xlab="Year",type='n',xaxs='i',axes=F)
	xx <- c(yr, rev(yr))
	yy <- c(sb[1,],rev(sb[3,]))
	polygon(xx, yy, col=shadecol[1], border='NA')
	lines(yr, sb[2,],lty=1,col=linecol[1],lwd=2)
	box()
	axis(side=1); axis(side=2)
	legend("topright", legend = "Median with 95% Interval", bty = 'n')
	legend("bottomleft", legend = paste("Median(Depletion) ", hist.yrs[length(hist.yrs)], " = ", round(sb[2,length(hist.yrs)],3)) , bty = 'n')
	dev.off()

	# OFL and ABC plots
	pngfun (file = "/OFL_ABC.png")
	par(mfrow=c(1,2),mar=c(4,4,4,4),oma=c(2,2,2,2))
	ymax = max(rep.list$OFL)
	boxplot(t(rep.list$OFL), ylim = c(0, ymax), ylab = "OFL", xlab ="Year")
	boxplot(t(rep.list$ForeCat), ylim = c(0, ymax), , ylab = "Adjusted Catch (ABC)", xlab ="Year")
	dev.off()

	# Plot the Q distribution
	for (i in 1:n.survey){
		pngfun (file = paste0("/Survey_Q_", i, ".png"))
		xx = length(val)
		n = matchfun(string = "Survey_Q", obj = colnames(val[[xx]])) + i - 1
		hist(exp(val[[xx]][,n]), xlim = c(0, max(exp(val[[xx]][,n]))), main = "", xlab = paste0("Survey Q ", i))
		abline(v = median(exp(val[[xx]][,n])), col = 'red', lwd =2)
		legend("topright", legend = paste("Median = ",round(median(exp(val[[xx]][,n])),3)) , bty = 'n')
		dev.off()
	}

	# Plot the estimate extra variance for the survey
	if(n.extra.se > 0 ){
	for(i in 1:n.extra.se){
		pngfun (file = paste0("/Added_Survey_Variance_",i,".png"))
		xx = length(val)
		n = matchfun(string = "extra_se", obj = colnames(val[[xx]])) + i - 1
		hist(val[[xx]][,n], xlim = c(0, max(val[[xx]][,n])), main = "", xlab = paste0("Added Survey Variance", i))
		abline(v = median(val[[xx]][,n]), col = 'red', lwd =2)
		legend("topright", legend = paste("Median = ",round(median(val[[xx]][,n]),3)) , bty = 'n')
		dev.off()
	}}


	# Write output tables
	sb.ci   = apply(rep.list$SB, 1, quantile, c(0.025, 0.975))
	tb.ci   = apply(rep.list$TotBio, 1, quantile, c(0.025, 0.975))
	smry.ci = apply(rep.list$SmryBio, 1, quantile, c(0.025, 0.975))
	d.ci    = apply(rep.list$Bratio, 1, quantile, c(0.025, 0.975))
	spr.ci  = apply(rep.list$SPR, 1, quantile, c(0.025, 0.975))
	exp.ci  = apply(rep.list$Exploitation, 1, quantile, c(0.025, 0.975))

	out     = cbind(comma(apply(rep.list$SB, 1, median), digits = 0),
				paste0(comma(sb.ci[1,],digits = 0), "\u2013", comma(sb.ci[2,],digits = 0)),

				comma(apply(rep.list$TotBio,1, median), digits = 0),
				paste0(comma(tb.ci[1,],digits = 0), "\u2013", comma(tb.ci[2,],digits = 0)),

				comma(apply(rep.list$SmryBio,1, median), digits = 0),
				paste0(comma(smry.ci[1,],digits = 0), "\u2013", comma(smry.ci[2,],digits = 0)),

				comma(apply(rep.list$Bratio, 1, median), digits = 2)),
				c("-", paste0(comma(d.ci[1,],digits = 2), "\u2013", comma(d.ci[2,],digits = 2))) 

				comma(apply(rep.list$SPR, 1, median), digits = 3)),
				c("-", paste0(comma(spr.ci[1,],digits = 2), "\u2013", comma(spr.ci[2,],digits = 2))) 

				comma(apply(rep.list$Explotation, 1, median), digits = 3)),
				c("-", paste0(comma(exp.ci[1,],digits = 3), "\u2013", comma(exp.ci[2,],digits = 3))) )

	colnames(out) = c("SB", "95%", "Total_Biomass", "95%", "Summary_Biomass", "95%","Depletion", "95%", "SPR", "95%", "Exploitation", "95%")
	write.csv(out, file = paste0(dir, "/Median_TimeSeries.csv"))

	ofl.ci = apply(rep.list$OFL, 1, quantile, c(0.025, 0.975))
	abc.ci = apply(rep.list$ABC, 1, quantile, c(0.025, 0.975))
	out = cbind(comma(apply(rep.list$OFL, 1, median), digits = 0),
				paste0(comma(ofl.ci[1,],digits = 0), "\u2013", comma(ofl.ci[2,],digits = 0)),
				comma(apply(rep.list$ABC, 1, median), digits = 0),
				paste0(comma(abc.ci[1,],digits = 0), "\u2013", comma(abc.ci[2,],digits = 0)) )
	colnames(out) = c("OFL", "95%", "ABC", "95%")
	write.csv(out, file = paste0(dir, "/Median_OFL_ABC.csv"))
}