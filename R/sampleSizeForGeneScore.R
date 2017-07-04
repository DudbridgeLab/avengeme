#' Sample size calculations for polygenic scores
#'
#' Calculates the size of training sample to achieve a given AUC, R2 or power in the target sample.
#'
#' The sample size is estimated by numerical optimisation.  For each possible sample size,
#' the P-value threshold is identified for selecting markers into the polygenic score, such that targetQuantity is maximised.
#'
#' @param targetQuantity Either "AUC", "R2" or "power" (case insensitive).
#' @param targetValue The value of the targetQuantity for which to calculate sample size.
#' @param n2 Target sample size.  Only relevant when targetQuantity is "power". By default set equal to the training sample size.
#' @inheritParams polygenescore

#' @return A list with the following elements:
#' \itemize{
#'	\item{\code{n} Required sample size for the training sample.  This is the total sample size: to obtain the number of cases, multiply by the sampling fraction.}
#' 	\item{\code{p} P-value threshold for selecting markers into the polygenic score, such that the target value is achieved with the minimum sample size.}
#' 	\item{\code{max} Maximum targetQuantity possible if the training sample size were increased to infinity (actually 1e10).}
#' }

#' @examples
#' # AUC= 0.75 in breast cancer.  See Table 4, row 4, column 3 in Dudbridge (2013).
#' sampleSizeForGeneScore("AUC",0.75,nsnp=100000,vg1=0.44/2,pi0=0.90,binary=TRUE,
#' prevalence=0.036,sampling=0.5)
#' # $n
#' # [1] 313981.4
#' #
#' # $p
#' # [1] 0.007500909
#' #
#' # $max
#' # [1] 0.788842
#' #
#' # Number of cases
#' sampleSizeForGeneScore("AUC",0.75,nsnp=100000,vg1=0.44/2,pi0=0.90,binary=TRUE,
#' prevalence=0.036,sampling=0.5)$n/2
#' # [1] 156990.7

#' @author Frank Dudbridge

#' @references
#' Dudbridge F (2013) Power and predictive accuracy of polygenic risk scores. PLoS Genet 9:e1003348
#' @export
sampleSizeForGeneScore=function(targetQuantity,
                                targetValue,
                                nsnp,
                                n2=NA,
                                vg1=0,
                                cov12=vg1,
                                pi0=0,
                                weighted=TRUE,
                                binary=FALSE,
                                prevalence=0.1,
                                sampling=prevalence,
                                lambdaS=NA,
                                shrinkage=FALSE,
                                logrisk=FALSE,
                                alpha=0.05,
					  r2gx=0,
					  corgx=0,
					  r2xy=0,
					  adjustedEffects=FALSE) {

###
# internal objective function, returns AUC/R2/Power for given pupper and n1
 obj2=function(p,n1) {
    if (is.na(n2)) {n2here=n1}
    else n2here=n2;
    pgs=polygenescore(nsnp,c(n1,n2here),vg1=vg1,cov12=cov12,pi0=pi0,pupper=c(0,p),weighted=weighted,binary=binary,prevalence=prevalence,sampling=sampling,
      lambdaS=lambdaS,shrinkage=shrinkage,logrisk=logrisk,alpha=alpha,r2gx=r2gx,corgx=corgx,r2xy=r2xy,adjustedEffects=adjustedEffects)
    if (tolower(targetQuantity)=="auc") {return(pgs$AUC)}
    if (tolower(targetQuantity)=="r2") {return(pgs$R2)}
    if (tolower(targetQuantity)=="power") {return(pgs$power)}
  }
### end of obj2

###
# internal objective function, maximises AUC/R2/Power over pupper, given n1
#   then returns squared difference from target value
  obj1=function(n1) {
    (optimise(obj2,c(0,1),n1,maximum=T)$objective-targetValue)^2
  }
### end of obj1

# get the order of magnitude for the sample size
  maxN=1
  while (optimise(obj2,c(0,1),maxN,maximum=T)$objective<targetValue) {
    maxN=maxN*10
  }

# minimise squared difference bewteen AUC/R2/Power and the target value
  fit=optimise(obj1,c(maxN/10,maxN))

  return(list(n=fit$minimum,
              p=optimise(obj2,c(0,1),fit$minimum,maximum=T)$maximum,
              max=optimise(obj2,c(0,1),1e10,maximum=T)$objective))
}
