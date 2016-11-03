#########
# avengeme.R: functions for analysis of polygenic scores
# Frank Dudbridge
# December 2013 - July 2015
# February 2014
#
# Additive Variance Explained and Number of Genetic Effects Method of Estimation
#
# Full documentation is in progress but the comments in this file should help users to get started
# Any questions to frank.dudbridge@lshtm.ac.uk
#
# This file contains three functions:
# polygenescore(): power, area under curve and correlation for estimated gene scores
# sampleSizeForGeneScore(): calculate the size of training sample required to achieve a given AUC, R2 or power
# estimatePolygenicModel(): estimate parameters of genetic model from a set of p-values obtained for the polygenic score in the target sample
#
# the functions have a common set of core parameters, which are described below under polygenescore
#
#
#########
# polygenescore: power, area under curve and correlation for estimated gene scores
#
#
# INPUTS
# nsnp = number of independent SNPs in the gene score
# n = vector containing training and target sample sizes
# vg1 = proportion of total variance that is explained by genetic effects in training sample
# cov12 = covariance between genetic effect sizes in the two populations
#         this must be <=sqrt(vg1)
#         if the genetic effects are fully correlated, set cov12=vg1
#         or leave unspecified, as this is the default
# pi0 = proportion of SNPs with no effects on training trait
# pupper = vector of thresholds on p-value for selection from training sample
#          first element is the lower bound of the first interval
#          second element is the upper bound of the first interval
#          third element is the upper bound of the second interval, etc
# nested = T if the p-value intervals are nested, ie have the same lower bound
#            if false, lower bound of the 2nd interval = upper bound of the 1st
#            and so on. Default nested=T
# weighted = T if effect sizes used as weights in forming gene score
#              if false, unweighted score is used, ie sum of risk alleles
#              Default nested = T
# binary = T if training trait is binary
#          by default, binary status of the target< trait is the same as
#          the training trait.  If different, make binary into a vector
#          with 2 elements          
# prevalence = disease prevalence in training sample
#              by default, prevalence is the same in the target sample
#              If different, make prevalence a vector with 2 elements
# sampling = case/control sampling fraction in training sample
#            by default, equals the prevalence, as in a cohort study
#            if sampling fraction is different in the target sample,
#            make sampling into a vector with 2 elements
# lambdaS = sibling relative recurrence risk in training sample, can be specified instead of vg1
# shrinkage = T if effect sizes are to be shrunk to BLUPs
# logrisk = T if binary trait arises from log-risk model
#           rather than liability threshold
# alpha = type-1 error for testing association in target sample

# OUTPUTS
# R2 = squared correlation between estimated gene score and target trait
# NCP = non-centrality parameter of association test between score and target trait
# p = expected p-value of association test
# power = power of association test
# AUC = for binary traits, area under ROC curve
# MSE = for quantitative traits, mean square error between target trait
#       and estimated gene score
# error = error message, if any
#########

polygenescore=function(nsnp,
                       n,
                       vg1=0,
                       cov12=vg1,
                       pi0=0,
                       pupper=c(0,1),
		       nested=T,
                       weighted=T,
                       binary=c(F,F),
                       prevalence=c(0.1,0.1),
                       sampling=prevalence,
                       lambdaS=NA,
                       shrinkage=F,
                       logrisk=F,
                       alpha=0.05) {

errorMsg=""
if (min(pupper)<0 | max(pupper)>1)
  errorMsg="Error: entries in pupper should be in (0,1)"
if (errorMsg!="")
  return(list(R2=NULL,NCP=NULL,p=NULL,power=NULL,AUC=NULL,MSE=NULL,error=errorMsg))

if (length(n)==1) n=rep(n,2)
if (length(binary)==1) binary=rep(binary,2)
if (length(prevalence)==1) prevalence=rep(prevalence,2)
if (length(sampling)==1) sampling=rep(sampling,2)

# variance of the phenotype
varY1 = 1
if (binary[1]) {
  varY1 = sampling[1]*(1-sampling[1])
}

varY2 = 1
if (binary[2]) {
  varY2 = sampling[2]*(1-sampling[2])
}

# sampling variance of each beta
samplingVar = varY1/n[1]

# conversion from lambdaS to vg
if (!is.na(lambdaS)) {
  if (logrisk) {
    vg1=2*log(lambdaS)
  }
  else {
    t=-qnorm(prevalence[1])
    t1=qnorm(1-lambdaS*prevalence[1])
    vg1 = 2*(t-t1*sqrt(1-(t^2-t1^2)*(1-t*prevalence[1]/dnorm(t))))/(dnorm(t)/prevalence[1]+t1^2*(dnorm(t)/prevalence[1]-t))
    if (vg1>1 | is.na(vg1)) vg1=1
  }
}

# variance of nonzero betas
betaVar = vg1/(nsnp*(1-pi0))
betaCov = cov12/(nsnp*(1-pi0))

# transform from liability scale to observed scale
if (logrisk) {
  liab2obs=prevalence*sampling*(1-sampling)/prevalence/(1-prevalence)
}
else {
  liab2obs=dnorm(qnorm(prevalence))*sampling*(1-sampling)/prevalence/(1-prevalence)
}

if (binary[1]) betaVar = betaVar*liab2obs[1]^2
if (binary[1]) betaCov = betaCov*liab2obs[1]
if (binary[2]) betaCov = betaCov*liab2obs[2]

shrink=1
if (shrinkage) {
  shrink = 1-samplingVar/(betaVar*(1-pi0)+samplingVar)
#  betaVar = betaVar*shrink^2
#  betaVar2 = betaVar2*shrink^2
#  samplingVar = samplingVar*shrink^2
}

R2=c()
NCP=c()
p=c()
power=c()
AUC=c()
MSE=c()

# threshold on betahat based on its p-value
  plower=c(pupper[-length(pupper)])
  pupper=pupper[-1]
  if (nested) plower=rep(plower[1],length(pupper))

  betaHatThreshLo = -qnorm(plower/2)*sqrt(samplingVar)
  betaHatThreshHi = -qnorm(pupper/2)*sqrt(samplingVar)

# expected number of selected SNPs
  betaHatSD = sqrt(betaVar+samplingVar)
  probTruncBeta = 2*nsnp*(1-pi0)*abs(pnorm(-betaHatThreshHi,sd=betaHatSD)-
                                     pnorm(-betaHatThreshLo,sd=betaHatSD))
  nullHatSD = sqrt(samplingVar)
  probTruncNull = 2*nsnp*pi0*abs(pnorm(-betaHatThreshHi,sd=nullHatSD)-
                                 pnorm(-betaHatThreshLo,sd=nullHatSD))

# variance of the estimated gene score
  if (weighted) {
    term1=rep(0,length(plower))
    term2=rep(0,length(plower))
    w=which(plower>0)
    term1[w]=betaHatThreshLo[w]/betaHatSD*dnorm(betaHatThreshLo[w]/betaHatSD)
    w=which(pupper>0)
    term2[w]=betaHatThreshHi[w]/betaHatSD*dnorm(betaHatThreshHi[w]/betaHatSD)
    varBetaHat = betaHatSD^2*(1+(term1-term2)/(pnorm(betaHatThreshHi/betaHatSD)-pnorm(betaHatThreshLo/betaHatSD)))
    w=which(plower>0)
    term1[w]=betaHatThreshLo[w]/nullHatSD*dnorm(betaHatThreshLo[w]/nullHatSD)
    w=which(pupper>0)
    term2[w]=betaHatThreshHi[w]/nullHatSD*dnorm(betaHatThreshHi[w]/nullHatSD)
    varNullHat = samplingVar*(1+(term1-term2)/(pnorm(betaHatThreshHi/nullHatSD)-pnorm(betaHatThreshLo/nullHatSD)))
    varGeneScoreHat = varBetaHat*probTruncBeta+varNullHat*probTruncNull
  }
  else {
    varGeneScoreHat = probTruncBeta+probTruncNull
  }

# covariance between Y2 and estimated gene score
  if (weighted) {
# regression coefficient in SNPs with effects on both traits
    scoreCovariance = betaCov/(betaVar+samplingVar)
# covariance in SNPs with effects
    scoreCovariance = scoreCovariance*varBetaHat*probTruncBeta
  }
  else {
    scoreCovariance=rep(0,length(plower))
    for(i in 1:length(plower)) {
      scoreCovariance[i] = 2*betaCov/betaVar*(1-pi0)*nsnp*
        integrate(discordantSign,0,Inf,sqrt(betaVar),betaHatThreshLo[i],betaHatThreshHi[i],sqrt(samplingVar),abs.tol=1e-12,stop.on.error=F)$value
    }
  }

# Coefficient of determination!
  R2 = scoreCovariance^2/varGeneScoreHat/varY2
# catch impossible configurations
  if (is.nan(max(R2)) | max(R2)>1 | vg1>0 & abs(cov12)/sqrt(vg1)>1) {
    mynan=rep(NaN,length(plower))
    return(list(R2=mynan,NCP=mynan,p=mynan,power=mynan,AUC=mynan,MSE=mynan,error="cov12 incompatible with vg1"))
  }

# Non-centrality parameter!
  NCP=n[2]*R2/(1-R2)
# Power!
  power=pchisq(qchisq(1-alpha,1),1,lower=F,ncp=NCP)

  thresholdDensity = dnorm(qnorm(prevalence[2]))/prevalence[2]
  caseMean = thresholdDensity*R2*varY2/liab2obs[2]^2
  caseVariance = R2*varY2/liab2obs[2]^2*(1-caseMean*(thresholdDensity+qnorm(prevalence[2])))
  thresholdDensity = dnorm(qnorm(prevalence[2]))/(1-prevalence[2])
  controlMean = -thresholdDensity*R2*varY2/liab2obs[2]^2
  controlVariance = R2*varY2/liab2obs[2]^2*(1+controlMean*(thresholdDensity-qnorm(prevalence[2])))

# debugging
#print(c(probTruncBeta,probTruncNull,varGeneScoreHat,scoreCovariance,caseMean,controlMean,caseVariance,controlVariance))
#print(varGeneScoreHat)

# area under ROC curve!
  if (binary[2]) {
    if (logrisk) {
      AUC = pnorm(sqrt(R2*(1-prevalence[2])^2/sampling[2]/(1-sampling[2])/2))
    }
    else {
      AUC = pnorm((caseMean-controlMean)/sqrt(caseVariance+controlVariance))
    }
    MSE=NULL
  }
  else {
    AUC = NULL
    MSE = 1+shrink^2*varGeneScoreHat-2*shrink*scoreCovariance
  }

# R2 on liability scale for binary traits
  if (binary[2]) R2 = R2/liab2obs[2]^2*sampling[2]*(1-sampling[2])

return(list(R2=R2,NCP=NCP,p=pchisq(NCP+1,1,lower=F),power=power,AUC=AUC,MSE=MSE,error=""))
}

# numerical integration in covariance for unweighted score
discordantSign=function(x,xsigma,threshLo,threshHi,asigma) {
  x*dnorm(x,sd=xsigma)*
  (pnorm(threshLo,mean=x,sd=asigma)-pnorm(threshHi,mean=x,sd=asigma)-pnorm(-threshHi,mean=x,sd=asigma)+pnorm(-threshLo,mean=x,sd=asigma))
  
}
#########
#end of polygenescore
#########

#########
# sampleSizeForGeneScore: calculate the size of training sample for a given AUC, R2 or power
# INPUTS:
# targetQuantity: "AUC", "R2" or "Power"
# targetValue: the value of AUC/R2/Power for which to calculate sample size
# n2: target sample size, by default set equal to the training sample size
# Other parameters as for polygenescore
#   but note that pupper is not used
# 
# OUTPUTS:
# n = required sample size for the training sample
# p = p-value threshold for selecting SNPs in the score, such that the target value is achieved with the minimum sample size
# max = maximum AUC/R2/Power possible if the training sample contained all living humans (n1=1e10)
#########
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
#########
#end of sampleSizeForGeneScore
#########


#########
# estimatePolygenicModel: estimate parameters of genetic model from a set of p-values for the polygenic score tested in the target sample
# each p-value corresponds to a polygenic score constructed from SNPs whose p-values in the training data are in a given interval
# possible parameters to estimate are vg1, cov12, vg2, pi01, pi02, in any combination
# a parameter will be estimated if its input value is unspecified or NA

# INPUTS:
# p = vector of p-values (or Z statistics) for polygenic scores tested in the target data
#   automatically detects Z statistics if some entries of p are greater than 1
# vg = proportion of total variance that is explained by genetic effects in training sample
#      by default, the variance is the same for the target sample
#      to specify a different variance, make vg a 2-element vector
# cov12 = covariance between genetic effect sizes in the two populations
# pi0 = proportion of SNPs with no effects on training trait
#       by default, the proportion is the same for the target trait
#       to specify a different proportion, make pi0 a 2-element vector
# boot = number of bootstrap replicates to estimate approximate confidence interval
#        if boot==0 (default), an analytic interval is calculated using profile likelihood
#        if boot>0, a bootstrap interval is given
#        these intervals assume that the input p-values are independent
#        this assumption is generally untrue and the interval will be slightly smaller than it should be
# bidirectional = T if results are also given when exchanging the training and target samples
#               in this case, vg2 can also be estimated, as can pi0 in the target sample
#               in this case, input vector p should be twice as long
#               with the list of p-values for training/target followed by
#               the list for target/training
# fixvg2pi02 = T if the same genetic model is assumed for the training and target samples
#                this fixes the target variance and the covariance to both equal the variance in the training sample
#                also fixes the proportion of null SNPs in the target sample to equal the proportion in the training sample
# other parameters as in polygenescore
#########
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
                                alpha=0.05,
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
    errorMsg=paste("Error: number of free parameters (",nparam,") is less than number of results (",length(p),")",sep="")


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
      bootNCP[1:(length(p)/2)]=polygenescore(nsnp,n,vghere[1],cov12here,pi0here[1],pupper,nested,weighted,binary,prevalence,sampling,lambdaS[1],shrinkage,logrisk,alpha)$NCP
      bootNCP[(length(p)/2+1):length(p)]=polygenescore(nsnp,rev(n),vghere[2],cov12here,pi0here[2],pupper,nested,weighted,rev(binary),rev(prevalence),rev(sampling),lambdaS[2],shrinkage,logrisk,alpha)$NCP
    }
    else {
      bootNCP=polygenescore(nsnp,n,vghere[1],cov12here,pi0here[1],pupper,nested,weighted,binary,prevalence,sampling,lambdaS[1],shrinkage,logrisk,alpha)$NCP
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

#########
#end of estimatePolygenicModel
#########
