#' Plot predictive accuracy
#'
#' \code{plotAccuracy} plots the AUC or the R2 as a function of training sample size.
#'
#' AUC is plotted for binary traits, R2 for quantitative traits.  At each point, the p-value threshold is identified for selecting markers into the polygenic score,
#' such that the AUC or R2 is maximised.

#' @param xlim Vector of 2 elements, giving the range of sample size to display on the x-axis, in 1000s.  For binary traits this is the number of cases.
#' @param ylim Range of AUC/R2 to display on y-axis.
#' @param nsnp Number of independent SNPs in the gene score.
#' @param fix TRUE if the same genetic model is assumed for the training and target samples.
#' @param plot TRUE is a new plot is to be drawn, otherwise draw lines on the existing plot.
#' @param col Colour in which to plot.
#' @param breakeven Value of AUC/R2 for which the minimum sample size will be estimated.
#' @param lty Line type parameter for R plots.
#' @inheritParams polygenescore

#' @import graphics

#' @return A list with the following elements:
#' \itemize{
#'	\item{\code{limit} Value of AUC/R2 at the maximum sample size plotted.}
#'	\item{\code{breakeven} Sample size at which the AUC/R2 exceeds the value specified by the breakeven parameter}
#'	\item{\code{plimit} Optimal p-value threshold at the maximum sample size plotted}
#' }

#' @examples
#' # Breast cancer with 90% null markers, from figure 3 in Dudbridge (2013)
#' plotAccuracy(vg1=0.44/2,pi0=0.90,fix=TRUE,binary=TRUE,prevalence=0.036)

#' @author Frank Dudbridge

#' @references
#' Dudbridge F (2013) Power and predictive accuracy of polygenic risk scores. PLoS Genet 9:e1003348

#' @export
plotAccuracy=function(xlim=c(1,500),ylim=0,nsnp=100000,vg1=0,pi0=0,cov12=NA,fix=T,binary=F,prevalence=0,sampling=.5,r2gx=0,corgx=0,r2xy=0,adjustedEffects=F,plot=T,col="black",breakeven=0.5,lty=1) {
  obj=function(p) {if (binary) n=1000/sampling else n=1000;a=polygenescore(n=n*i,nsnp=nsnp,vg1=vg1,cov12=cov12,prevalence=prevalence,binary=binary,sampling=sampling,pi0=pi0,pupper=c(0,p),r2xy=r2xy,r2gx=r2gx,corgx=corgx,adjustedEffects=adjustedEffects,riskthresh=0.0);if (binary) 1-a$AUC else 1-a$R2}
  if (is.na(cov12)) cov12=vg1
  p1=NULL
  p2=NULL
  x=seq(xlim[1],xlim[2],1)
  for(i in x) {
    fit=optimise(obj,c(0,1))
    p1=c(p1,1-fit$obj)
    p2=c(p2,fit$min)
  }
  if (plot)
    if (binary) {
      if (length(ylim)!=2) ylim=c(0.5,1)
      plot(x,p1,type="l",ylim=ylim,xlim=xlim,ylab="AUC",xlab="n (1000 cases)",lty=lty)
    }
    else {
      if (length(ylim)!=2) ylim=c(0,1)
      plot(x,p1,type="l",ylim=ylim,xlim=xlim,ylab="R2",xlab="n (1000 subjects)",lty=lty)
    }
  else
    lines(x,p1,col=col,lty=lty)
  return(list(limit=p1[length(p1)],breakeven=x[which(p1>=breakeven)[1]],plimit=p2[length(p2)]))
}

#########
#end of plotAccuracy
#########
