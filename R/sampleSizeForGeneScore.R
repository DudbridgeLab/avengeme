#' Sample Size for Genetic Score
#' Calculate the size of training sample for a given AUC, R2 or power targetQuantity: "AUC", "R2" or "Power"
#' @param targetValue: the value of AUC/R2/Power for which to calculate sample size
#' @param n2: target sample size, by default set equal to the training sample size
#' @return 
#' \itemize{
#'	\item{"n = required sample size for the training sample"}
#' 	\item{"p = p-value threshold for selecting SNPs in the score, such that the target value is achieved with the minimum sample size"}
#' 	\item{"max = maximum AUC/R2/Power possible if the training sample contained all living humans (n1=1e10)"}
#' }
sampleSizeForGeneScore=function(targetQuantity,
                                targetValue,
                                nsnp,
                                n2=NA,
                                vg1=0,
                                cov12=vg1,
                                pi0=0,
                                weighted=T,
                                binary=F,
                                prevalence=0.1,
                                sampling=prevalence,
                                lambdaS=NA,
                                shrinkage=F,
                                logrisk=F,
                                alpha=0.05) {

###
# internal objective function, returns AUC/R2/Power for given pupper and n1 
 obj2=function(p,n1) {
    if (is.na(n2)) {n2here=n1}
    else n2here=n2;
    pgs=polygenescore(nsnp,n1,vg1,n2here,cov12,pi0,plower=0,pupper=p,weighted,binary,prevalence,sampling,lambdaS,shrinkage,logrisk,alpha)
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
