#' Polygenic Score:
#' Power, area under curve and correlation for estimated gene scores.
#'
#' @param nsnp number of independent SNPs in the gene score
#' @param n vector containing training and target sample sizes
#' @param vg1 proportion of total variance that is explained by genetic effects in training sample
#' @param cov12 covariance between genetic effect sizes in the two populations. This must be <=sqrt(vg1) if the genetic effects are fully correlated, set cov12=vg1 or leave unspecified, as this is the default.
#' @param pi0 proportion of SNPs with no effects on training trait
#' @param pupper vector of thresholds on p-value for selection from training sample. first element is the lower bound of the first interval, second element is the upper bound of the first interval third element is the upper bound of the second interval, etc
#' @param nested T if the p-value intervals are nested, ie have the same lower bound if false, lower bound of the 2nd interval = upper bound of the 1st and so on. Default nested=T
#' @param weighted T if effect sizes used as weights in forming gene score if false, unweighted score is used, ie sum of risk alleles Default nested = T
#' @param binary T if training trait is binary by default, binary status of the target< trait is the same as the training trait.  If different, make binary into a vector with 2 elements          
#' @param prevalence disease prevalence in training sample by default, prevalence is the same in the target sample If different, make prevalence a vector with 2 elements
#' @param sampling case/control sampling fraction in training sample by default, equals the prevalence, as in a cohort study if sampling fraction is different in the target sample, make sampling into a vector with 2 elements
#' @param lambdaS sibling relative recurrence risk in training sample, can be specified instead of vg1
#' @param shrinkage T if effect sizes are to be shrunk to BLUPs
#' @param logrisk T if binary trait arises from log-risk model rather than liability threshold
#' @param alpha type-1 error for testing association in target sample

#' @return
#' \itemize{
#'	\item{"R2 = squared correlation between estimated gene score and target trait"}
#'	\item{"NCP = non-centrality parameter of association test between score and target trait"}
#'	\item{"p = expected p-value of association test"}
#'	\item{"power = power of association test"}
#'	\item{"AUC = for binary traits, area under ROC curve"}
#'	\item{"MSE = for quantitative traits, mean square error between target trait and estimated gene score"}
#'	\item{"error = error message, if any"}
#' }

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


