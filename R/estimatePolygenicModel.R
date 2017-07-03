#' Estimate polygenic model
#'
#' \code{estimatePolygenicModel} estimates the parameters of an underlying genetic model from the results of association tests of a polygenic score.
#'
#' The input is a vector of P-values or (signed) Z-statistics from the association test of the polygenic score in the target sample.
#' P-values are assumed if all the values are in (0,1), otherwise Z-statistics are assumed.
#' Each P-value corresponds to the association test of a polygenic score consisting of SNPs with training sample P-values in a specific interval.
#' Up to five parameters can be estimated: vg[1], cov12, vg[2], pi0[1], pi0[2].
#' The number of input P-values must be greater than or equal to the number of estimated parameters, otherwise an error message is returned.
#' Any combination of parameters can be estimated.  A parameter will be estimated if its input value is unspecified or \code{NA}.

#' @param p Vector of P-values or Z-statistics for polygenic scores tested in the target data. Automatically detects Z-statistics if some entries of p are greater than 1 or less than 0.
#' @param vg Proportion of variance explained by genetic effects in training sample. By default, the variance explained is the same in the target sample; otherwise vg should be a vector with two elements for the training and target samples respectively.
#' @param pi0 Proportion of markers with no effect on the training trait.  By default, the proportion is the same for the target trait; otherwise pi0 should be a vector with two elements for the training and target samples respectively.
#' @param boot Number of bootstrap replicates to estimate approximate confidence intervals. If boot==0 (default), an analytic interval is calculated using profile likelihood.
#' if boot>0, a bootstrap interval is estimated. These intervals assume that the input P-values are independent; this assumption is generally untrue and the interval will be slightly smaller than it should be.
#' @param bidirectional TRUE if results are also given when exchanging the role of training and target samples.
#' In this case, vg2 can also be estimated, as can pi0 in the target sample. The input vector p should now be twice as long with the list of P-values for training/target followed by the list for target/training.
#' @param initial Specify starting values for numerical maximisation of the likelihood.  The number of elements must equal the number of estimated parameters, and follows the order vg[1], vg[2], pi0[1], pi0[2], cov12, for those parameters that are actually being estimated.  Default 0.5 for all parameters.
#' @param fixvg2pi02 TRUE if the same genetic model is assumed for the training and target samples.
#' This fixes the target variance and the covariance to both equal the variance explained in the training sample, vg1. Also fixes the proportion of null markers in the target sample to equal the proportion in the training sample.
#' @param option Parameter used in method development.  Default 0, fits the model by maximum likelihood for Z statistics.  1 and 2 fit the model by least squares to chisq and Z statistics respectively.  3 fits by maximum likelihood for chisq statistics.
#' @inheritParams polygenescore
#'
#' @return A list with elements corresponding to the estimated genetic model.
#' Values fixed at input are returned unchanged with a degenerate confidence interval.
#' Each element is a vector consisting of the point estimate followed by its lower and upper 95\% confidence limit.
#' \itemize{
#'	\item{\code{vg} Variance explained in the target trait.  If bidirectional estimation is selected, \code{vg} is a matrix with two rows corresponding to the training and target samples respectively.}
#'	\item{\code{cov12} Covariance between genetic effects in the two samples.}
#'	\item{\code{pi0} Proportion of markers with no effect on the training trait. If bidirectional estimation is selected, \code{pi0} is a matrix with two rows corresponding to the training and target samples respectively.}
#'	\item{\code{logLikelihood} Maximised log-likelihood at the fitted model.}
#'	\item{\code{error} Error message, if any}
#' }
#'
#' @examples
#' # Schizophrenia PGC2 study, rightmost column of table 5 in Palla & Dudbridge (2015)
#' # P-values from supplementary table 6 of Schizophrenia Working Group (2014)
#' # Other parameters as in table 1 of Palla & Dudbridge (2015)
#' pupper=c(0,5e-8,1e-6,1e-4,1e-3,0.01,0.05,0.1,0.2,0.5,1)
#' p=c(9.85087e-24, 4.44037e-36, 2.08048e-71, 8.0594e-103, 2.0587e-138,
#'     1.4131e-164,5.8954e-166,3.75e-164,7.9488e-159,2.3286e-157)
#' estimatePolygenicModel(p,103125,c(77195,5120),pupper=pupper,nested=TRUE,binary=TRUE,prevalence=0.01,sampling=c(0.425,0.515),fixvg2pi02=TRUE)
#' # $vg
#' # [1] 0.2449328 0.2352049 0.2547500
#' #
#' # $cov12
#' # [1] 0.2449328 0.2352049 0.2547500
#' #
#' # $pi0
#' # [1] 0.8520102 0.8354152 0.8669681
#' #
#' # $logLikelihood
#' # [1] 30.1139
#' #
#' # $error
#' # [1] ""
#'
#' # Genetic covariance between bipolar disorder and schizophrenia, table 6 in Palla & Dudbridge (2015)
#' # Nagelkerke R2 SCZ-BPD from table S5 of Cross Disorder Group (2013)
#' R2N=c(.0044,.0065,.015,.023,.024,.025,.024,.024,.025,.025)
#' n1=11922
#' p1=6664/n1 # sampling fraction
#' # Convert to observed scale R2
#' R2O=R2N*(1-p1^(2*p1)*(1-p1)^(2*(1-p1)))
#' # Convert to chisq statistics
#' X2=n*R2O/(1-R2O)
#' # Now the same for BPD-SCZ
#' R2N=c(0.002,0.0048,0.012,0.017,0.021,0.021,0.022,0.021,0.021,0.021)
#' n2=17012
#' p2=9032/n2
#' R2O=R2N*(1-p2^(2*p2)*(1-p2)^(2*(1-p2)))
#' X2=c(n2*R2O/(1-R2O),X2)
#' # Perform bidirectional estimation with Z-scores as the first argument
#' # Small difference from published result due to minor bug fixes
#' pupper=c(0,.0001,.001,.01,.05,.1,.2,.3,.4,.5,1)
#' estimatePolygenicModel(sqrt(X2),nsnp=83884,n=c(n1,n2),binary=TRUE,pupper=pupper,prevalence=c(0.01,0.01),sampling=c(p1,p2),bidirectional=TRUE)
#' # $vg
#' # vg      vgLo      vgHi
#' # [1,] 0.8800473 0.1784087 0.9999528
#' # [2,] 0.6339373 0.1258082 0.9999382
#' #
#' # $cov12
#' # [1] 0.2092269 0.1895850 0.2226666
#' #
#' # $pi0
#' # pi0     pi0Lo     pi0Hi
#' # [1,] 0.7297216 0.6453690 0.9138123
#' # [2,] 0.7237033 0.5936525 0.9127537
#' #
#' # $logLikelihood
#' # [1] 19.55162
#' #
#' # $error
#' # [1] ""
#'
#' @author Frank Dudbridge

#' @references
#' Palla L and Dudbridge F (2015) A fast method using polygenic scores to estimate the variance explained by genome-wide marker panels and the proportion of variants affecting a trait. Am J Hum Genet 97:250-259
#' @references
#' Schizophrenia Working Group of the Psychiatric Genomics Consoritum (2014) Biological insights from 108 schizophrenia-associated genetic loci.  Nature 511:421-427
#' @references
#' Cross Disorder Group of the Psychiatric Genomics Consortium (2013) dentification of risk loci with shared effects on five major psychiatric disorders: a genome-wide analysis.  Lancet 381:1371-1379
#' @export
estimatePolygenicModel=function(p,
                                nsnp,
                                n,
                                vg=c(NA,NA),
                                cov12=NA,
                                pi0=c(NA,NA),
                                pupper=1,
					  nested=T,
                                weighted=T,
                                binary=c(F,F),
                                prevalence=c(0.1,0.1),
                                sampling=prevalence,
                                lambdaS=c(NA,NA),
                                shrinkage=F,
                                logrisk=F,
                                option=0,
                                boot=0,
                                bidirectional=F,
                                initial=c(),
                                fixvg2pi02=F
                                ) {

# inverse logit function
  expit=function(x) {
    y=x
    w=which(abs(x)<500)
    y[w]=exp(x[w])/(1+exp(x[w]))
    w=which(x>500)
    y[w]=1
    w=which(x<(-500))
    y[w]=0
    y
  }

###
# internal objective function, returns the log-likelihood of the z-scores for the proposed model
# vg1here is a local copy of vg1 derived from the input parameter
  obj1=function(param,vg,cov12,pi0,returnCov12=F) {
    vghere=vg
    pi0here=pi0
    if (length(param)>1) param=expit(param)
    i=1
    if (is.na(vg[1])) {vghere[1]=param[i];i=i+1}
    if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {vghere[2]=param[i];i=i+1}
    if (is.na(pi0[1])) {pi0here[1]=param[i];i=i+1}
    if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {pi0here=c(pi0here[1],param[i]);i=i+1}
    if (fixvg2pi02) {
      vghere=rep(vghere[1],2)
      pi0here=rep(pi0here[1],2)
    }

    upperBound=sqrt(vghere[1])
    if (bidirectional) upperBound=min(upperBound,sqrt(vghere[2]))
    if (bidirectional) upperBound=min(upperBound,sqrt(vghere[1]*vghere[2])*(1-max(pi0here))/sqrt((1-pi0here[1])*(1-pi0here[2])))

    optCov12=function(cov12here) {
      intervalNCP=polygenescore(nsnp,n,vghere[1],cov12here,pi0here[1],pupper,nested,weighted,binary,prevalence,sampling,lambdaS[1],shrinkage,logrisk)$NCP
      if (bidirectional) intervalNCP=c(intervalNCP,polygenescore(nsnp,rev(n),vghere[2],cov12here,pi0here[2],pupper,nested,weighted,rev(binary),rev(prevalence),rev(sampling),lambdaS[2],shrinkage,logrisk)$NCP)

      if (option==0) {
        intervalNCP=sqrt(intervalNCP)
        if (cov12here<0) intervalNCP=-intervalNCP
      }

# other ways we considered for fitting the model
      if (option==3) s=-sum(dchisq(ncp^2,1,ncp=intervalNCP,log=T),na.rm=T)
      if (option==1) s=sum((intervalNCP-ncp^2)^2)
      if (option==2) s=sum((sqrt(intervalNCP)-ncp)^2)
      if (option==0) s=-sum(dnorm(ncp,mean=intervalNCP,log=T),na.rm=T)
      if (is.nan(s) | is.nan(max(intervalNCP))) s=1e10
      s
    }

    if (fixvg2pi02) result=optCov12(vghere[1])
    else {
      if (is.na(cov12)) {
        if (returnCov12) {
          if (option==0) result=optimise(optCov12,c(-upperBound,upperBound))$minimum
          else result=optimise(optCov12,c(0,upperBound))$minimum
        }
        else {
          if (option==0) result=optimise(optCov12,c(-upperBound,upperBound))$objective
          else result=optimise(optCov12,c(-upperBound,upperBound))$objective
        }
      }
      else result=optCov12(cov12)
    }
    result
  }
### end of obj1

### start of main function
if (length(n)==1) n=rep(n,2)
if (length(vg)==1) vg=rep(vg,2)
if (length(lambdaS)==1) lambdaS=rep(lambdaS,2)
if (length(pi0)==1) pi0=rep(pi0,2)
if (length(binary)==1) binary=rep(binary,2)
if (length(prevalence)==1) prevalence=rep(prevalence,2)
if (length(sampling)==1) sampling=rep(sampling,2)

# error messages
  errorMsg=""
  if (bidirectional & length(p)!=(length(pupper)-1)*2 | !bidirectional & length(p)!=length(pupper)-1)
    errorMsg=paste("Error: number of results (",length(p),") does not match number of p-value intervals (",length(pupper)-1,")",sep="")
  if (min(pupper)<0 | max(pupper)>1)
   errorMsg="Error: entries in pupper should be in (0,1)"
# number of parameters to be estimated
  nparam=is.na(vg[1])+is.na(pi0[1])
  if (bidirectional & !fixvg2pi02) nparam=nparam+is.na(vg[2])+is.na(pi0[2])
  if (!fixvg2pi02) nparam=nparam+is.na(cov12)
  if (nparam>length(p))
    errorMsg=paste("Error: number of free parameters (",nparam,") is greater than number of results (",length(p),")",sep="")


# convert input to z-statistics
  if (bidirectional) w=!is.nan(as.numeric(p[1:(length(p)/2)]) & !is.nan(as.numeric(p[(length(p)/2+1):length(p)])))
  else w=!is.nan(as.numeric(p))
  w=which(w)
  p=as.numeric(p[w])
  pupper=c(pupper[1],pupper[w+1])

  ncp=p
  if (max(p,na.rm=T)>1 | min(p,na.rm=T)<0) {
# trying out signed z-scores to allow negative correlation
    if (option!=0)    ncp=sqrt(p)
    if (option==0)    ncp=p
  }
  else {
    ncp=qnorm(p/2,lower=F)
  }
  ncp=as.numeric(ncp)

# default confidence intervals
  vgLo=vg
  vgHi=vg
  vgLo=vg
  vgHi=vg
  cov12Lo=cov12
  cov12Hi=cov12
  pi0Lo=pi0
  pi0Hi=pi0

  llhd=NULL

if (errorMsg=="" & nparam==0) {
  llhd=obj1(param=NULL,vg=vg,cov12=cov12,pi0=pi0)
}

if (errorMsg=="" & nparam>0) {

# estimate the free parameters
  if (nparam==1) {
    solution=optimise(obj1,c(0,1),vg=vg,cov12=cov12,pi0=pi0)
    llhd=solution$objective
    solution=solution$minimum
  }
  else {
    if (length(initial)>0) initial=log(initial/(1-initial))
    if (length(initial)<nparam) initial=c(initial,rep(0,nparam-length(initial)))
    if (length(initial)>nparam) initial=initial[1:nparam]
    solution=optim(initial,obj1,vg=vg,cov12=cov12,pi0=pi0)
    llhd=solution$value
    solution=expit(solution$par)
  }

# confidence intervals
  solutionLo=solution
  solutionHi=solution
  if (boot==0) {
    if (option==0 | option==3) {
# profile likelihood confidence interval
    if (nparam==1) {
      if (!fixvg2pi02 & is.na(cov12)) {
        solution=obj1(0,vg,NA,pi0,T)
# objective function to find the parameter whose deviance is 3.841 from the maximum
        obj2=function(param) {
          (2*(obj1(param,vg=vg,cov12=param,pi0=pi0)-llhd)-qchisq(.95,1))^2
        }
        upperBound=sqrt(vg[1])
        if (bidirectional) upperBound=min(upperBound,sqrt(vg[1]*vg[2]))
        solutionLo=optimise(obj2,c(-upperBound,solution))$minimum
        solutionHi=optimise(obj2,c(solution,upperBound))$minimum
      }
      else {
# objective function to find the parameter whose deviance is 3.841 from the maximum
        obj2=function(param) {
          (2*(obj1(param,vg=vg,cov12=cov12,pi0=pi0)-llhd)-qchisq(.95,1))^2
        }
        solutionLo=optimise(obj2,c(0,solution))$minimum
        solutionHi=optimise(obj2,c(solution,1))$minimum
      }
    }
    else {
# objective function to find the parameter whose deviance is 3.841 from the maximum
# for each input parameter, the likelihood is maximised over the other nuisance parameters
        initial_tmp=initial
        obj2=function(param,name) {
          if (name=="vg1") vg[1]=param
          if (name=="vg2") vg[2]=param
          if (name=="pi01") pi0[1]=param
          if (name=="pi02") pi0[2]=param
          if (name=="cov12") cov12=param
          if (nparam==1) obj2llhd=optimise(obj1,interval=c(0,1),vg=vg,cov12=cov12,pi0=pi0)$objective
          else obj2llhd=optim(initial_tmp,obj1,vg=vg,cov12=cov12,pi0=pi0)$value
          (2*(obj2llhd-llhd)-qchisq(.95,1))^2
        }
      i=1
# profile CI for vg1
      if (is.na(vg[1])) {
        nparam=nparam-1
        initial_tmp=initial[-i]
        if (solution[i]>0) solutionLo[i]=optimise(obj2,c(0,solution[i]),name="vg1")$minimum else solutionLo[i]=0
        if (solution[i]<1) solutionHi[i]=optimise(obj2,c(solution[i],1),name="vg1")$minimum else solutionHi[i]=1
        i=i+1
        nparam=nparam+1
      }
# profile CI for vg2
      if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {
        nparam=nparam-1
        initial_tmp=initial[-i]
        if (solution[i]>0) solutionLo[i]=optimise(obj2,c(0,solution[i]),name="vg2")$minimum else solutionLo[i]=0
        if (solution[i]<1) solutionHi[i]=optimise(obj2,c(solution[i],1),name="vg2")$minimum else solutionHi[i]=1
        i=i+1
        nparam=nparam+1
      }
# profile CI for pi0[1]
      if (is.na(pi0[1])) {
        nparam=nparam-1
        initial_tmp=initial[-i]
        if (solution[i]>0) solutionLo[i]=optimise(obj2,c(0,solution[i]),name="pi01")$minimum else solutionLo[i]=0
        if (solution[i]<1) solutionHi[i]=optimise(obj2,c(solution[i],1),name="pi01")$minimum else solutionHi[i]=1
        i=i+1
        nparam=nparam+1
      }
# profile CI for pi0[2]
      if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {
        nparam=nparam-1
        initial_tmp=initial[-i]
        if (solution[i]>0) solutionLo[i]=optimise(obj2,c(0,solution[i]),name="pi02")$minimum else solutionLo[i]=0
        if (solution[i]<1) solutionHi[i]=optimise(obj2,c(solution[i],1),name="pi02")$minimum else solutionHi[i]=1
        i=i+1
        nparam=nparam+1
      }
# profile CI for cov12
      if (!fixvg2pi02 & is.na(cov12)) {
        nparam=nparam-1
        initial_tmp=initial[-i]
        vghere=vg
        pi0here=pi0
        j=1
        if (is.na(vg[1])) {vghere[1]=solution[j];j=j+1}
        if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {vghere[2]=solution[j];j=j+1}
        if (is.na(pi0[1])) {pi0here[1]=solution[j];j=j+1}
        if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {pi0here[2]=solution[j];j=j+1}
        cov12here=obj1(0,vghere,NA,pi0here,T)
        upperBound=1
        if (!bidirectional & !is.na(vg[1])) upperBound=sqrt(vghere[1])
        if (bidirectional) {
          if (!is.na(vg[1]) & is.na(vg[2])) upperBound=sqrt(vghere[1])
          if (is.na(vg[1]) & !is.na(vg[2])) upperBound=sqrt(vghere[2])
          if (!is.na(vg[1]) & !is.na(vg[2])) upperBound=sqrt(vghere[1]*vghere[2])
        }
        if (abs(cov12here)>0) {
          if (option==0) solutionLo[i]=optimise(obj2,c(-upperBound,cov12here),name="cov12")$minimum
          else solutionLo[i]=optimise(obj2,c(0,cov12here),name="cov12")$minimum
        } else solutionLo[i]=0
        if (abs(cov12here)<1) solutionHi[i]=optimise(obj2,c(cov12here,upperBound),name="cov12")$minimum else solutionHi[i]=cov12here/abs(cov12here)
        i=i+1
        nparam=nparam+1
      }
    }
    }
  }
  else {
# bootstrap confidence interval
    i=1
    vghere=vg
    pi0here=pi0
    cov12here=cov12
    if (is.na(vg[1])) {
      vghere[1]=solution[i]
      i=i+1
    }
    if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {
      vghere[2]=solution[i]
      i=i+1
    }
    if (is.na(pi0[1])) {
      pi0here[1]=solution[i]
      i=i+1
    }
    if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {
      pi0here[2]=solution[i]
      i=i+1
    }
    if (!fixvg2pi02 & is.na(cov12)) {
      cov12here=obj1(0,vghere,NA,pi0here,T)
      i=i+1
    }
    if (fixvg2pi02) {
      vghere=rep(vghere[1],2)
      cov12here=vghere[1]
      pi0here=rep(pi0here[1],2)
    }

    bootNCP=rep(0,length(p))
    if (bidirectional) {
      bootNCP[1:(length(p)/2)]=polygenescore(nsnp,n,vghere[1],cov12here,pi0here[1],pupper,nested,weighted,binary,prevalence,sampling,lambdaS[1],shrinkage,logrisk)$NCP
      bootNCP[(length(p)/2+1):length(p)]=polygenescore(nsnp,rev(n),vghere[2],cov12here,pi0here[2],pupper,nested,weighted,rev(binary),rev(prevalence),rev(sampling),lambdaS[2],shrinkage,logrisk)$NCP
    }
    else {
      bootNCP=polygenescore(nsnp,n,vghere[1],cov12here,pi0here[1],pupper,nested,weighted,binary,prevalence,sampling,lambdaS[1],shrinkage,logrisk)$NCP
    }
#    bootchisq=matrix(rchisq(boot*length(p),1,ncp=bootNCP),nrow=length(p))
    bootnorm=matrix(rnorm(boot*length(p),mean=sqrt(bootNCP)),nrow=length(p))
    if (cov12here<0) bootnorm=-bootnorm
    bootsolution=matrix(nrow=boot,ncol=nparam)
    ncp_orig=ncp
    for(i in 1:boot) {
#      ncp=sqrt(bootchisq[,i])
      ncp=bootnorm[,i]
      if (nparam==1) {
        if (!fixvg2pi02 & is.na(cov12)) bootsolution[i,1]=obj1(0,vg,NA,pi0,T)
        else bootsolution[i,1]=optimise(obj1,c(0,1),vg=vg,cov12=cov12,pi0=pi0)$minimum
      }
      else {
        bootsolution[i,]=expit(optim(initial,obj1,vg=vg,cov12=cov12,pi0=pi0)$par)
        fit=expit(optim(initial,obj1,vg=vg,cov12=cov12,pi0=pi0)$par)
        j=1
        if (is.na(vg[1])) {vghere[1]=fit[j];j=j+1}
        if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {vghere[2]=fit[j];j=j+1}
        if (is.na(pi0[1])) {pi0here[1]=fit[j];j=j+1}
        if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {pi0here[2]=fit[j];j=j+1}
        if (fixvg2pi02) {
          vghere=rep(vghere[1],2)
          cov12here=vghere[1]
          pi0here=rep(pi0here[1],2)
        }
        if (!fixvg2pi02 & is.na(cov12)) fit[j]=obj1(0,vghere,NA,pi0here,T)
        bootsolution[i,]=fit
      }
    }
    ncp=ncp_orig
    for(i in 1:nparam) {
      solutionLo[i]=quantile(bootsolution[,i],.025)
      solutionHi[i]=quantile(bootsolution[,i],.975)
    }
  }

# copy out final solutions
  i=1
  if (is.na(vg[1])) {
    vg[1]=solution[i]
    vgLo[1]=solutionLo[i]
    vgHi[1]=solutionHi[i]
    i=i+1
  }
  if (bidirectional & !fixvg2pi02 & is.na(vg[2])) {
    vg[2]=solution[i]
    vgLo[2]=solutionLo[i]
    vgHi[2]=solutionHi[i]
    i=i+1
  }
  if (is.na(pi0[1])) {
    pi0[1]=solution[i]
    pi0Lo[1]=solutionLo[i]
    pi0Hi[1]=solutionHi[i]
    i=i+1
  }
  if (bidirectional & !fixvg2pi02 & is.na(pi0[2])) {
    pi0[2]=solution[i]
    pi0Lo[2]=solutionLo[i]
    pi0Hi[2]=solutionHi[i]
    i=i+1
  }

  if (!fixvg2pi02 & is.na(cov12)) {
    cov12=obj1(0,vg,NA,pi0,T)
    cov12Lo=solutionLo[i]
    cov12Hi=solutionHi[i]
  }

  if (nparam>1) initial=expit(initial)
} # if errorMsg==""

  if (fixvg2pi02) {
    vg=rep(vg[1],2)
    vgLo=rep(vgLo[1],2)
    vgHi=rep(vgHi[1],2)
    cov12=vg[1]
    cov12Lo=vgLo[1]
    cov12Hi=vgHi[1]
    pi0=rep(pi0[1],2)
    pi0Lo=rep(pi0Lo[1],2)
    pi0Hi=rep(pi0Hi[1],2)
  }

 if (bidirectional) list(vg=cbind(vg,vgLo,vgHi),cov12=c(cov12,cov12Lo,cov12Hi),pi0=cbind(pi0,pi0Lo,pi0Hi),logLikelihood=llhd,error=errorMsg)
 else list(vg=c(vg[1],vgLo[1],vgHi[1]),cov12=c(cov12,cov12Lo,cov12Hi),pi0=c(pi0[1],pi0Lo[1],pi0Hi[1]),logLikelihood=llhd,error=errorMsg)

}


