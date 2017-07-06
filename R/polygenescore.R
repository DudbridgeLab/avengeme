#' Calculate power and predictive accuracy of a polygenic score
#'
#' Calculates measures of association for a polygenic score derived from a training sample
#' to predict traits in a target sample.
#'
#' The following setup is assumed. Two independent samples of genotypes are available; this could be
#' one sample of data split into two subsets.  One sample is termed the training sample, the other the
#' target sample.  Traits are measured in each sample; different traits could be measured in training and target
#' samples.  Subjects are assumed to be unrelated, and genotypes assumed to be
#' independent.  In practice we recommend LD-clumping methods, such as the --clump option in
#' PLINK, to ensure weak dependence between markers; in this case the methods are almost unbiased if an r2
#' threshold of 0.1 is used. Markers with P-values within a fixed range are selected from the training sample, and
#' then used to construct a polygenic score for each subject in the target sample.  The score can be tested for
#' association to the target trait, or used to predict individual trait values in the target sample.

#' @param nsnp Number of independent markers in the polygenic score.
#' @param n Vector with two elements, giving the total sizes of the training and target samples. In
#' case/control studies, n is the sum of the number of cases and number of controls. If only
#' one element of n is given, the training and target samples are assumed to be the same size.
#' No default - a value must be given
#' @param vg1 Proportion of variance explained by genetic effects in the training sample.
#' @param cov12 Covariance between genetic effect sizes in the two samples. If the effects are fully correlated then cov12<=sqrt(vg1).  If the effects are identical then cov12=vg1 (default).
#' @param pi0 Proportion of markers with no effect on the training trait.
#' @param pupper Vector of p-value thresholds for selecting markers from training sample. First element is the lower bound of the first interval, second element is the upper bound of the first interval, third element is the upper bound of the second interval, etc.
#' @param nested TRUE if the p-value intervals are nested, that is they have the same lower bound, which is the first element of pupper. If false, lower bound of the second interval is the upper bound of the first and so on.
#' @param weighted TRUE if estimated effect sizes are used as weights in forming the polygenic score. If false, an unweighted score is used, which is the sum of risk alleles carried.
#' @param binary TRUE if the training trait is binary. By default, the target trait is binary if the training trait is; otherwise binary should be a vector with two elements for the training and target samples respectively.
#' @param prevalence For a binary trait, prevalence in the training sample. By default, prevalence is the same in the target sample. Otherwise, prevalence should be a vector with two elements for the training and target samples respectively.
#' @param sampling For a binary trait, case/control sampling fraction in the training sample. By default, sampling equals the prevalence, as in a cohort study.  If the sampling fraction is different in the target sample, sampling should be a vector with two elements for the training and target samples respectively.
#' @param lambdaS Sibling relative recurrence risk in training sample, can be specified instead of vg1.
#' @param shrinkage TRUE if effect sizes are to be shrunk to BLUPs.
#' @param logrisk TRUE if binary trait arises from log-risk model rather than liability threshold.
#' @param alpha Significance level for testing association of the polygenic score in the target sample.
#' @param r2gx Proportion of variance in environmental risk score explained by genetic effects in training sample.
#' @param corgx Genetic correlation between environmental risk score and training trait.
#' @param r2xy Proportion of variance in training trait explained by environmental risk score.
#' @param adjustedEffects TRUE if polygenic and envrionmental scores are combined as a weighted sum. If FALSE, the scores are combined as an unweighted sum even if they are correlated.
#' @param riskthresh Absolute risk threshold for calculating net reclassification index.

#' @import stats

#' @return A list with elements containing quantities describing the association of the polygenic score with the target trait:
#' \itemize{
#'	\item{\code{R2} Squared correlation between polygenic score and target trait.}
#'	\item{\code{NCP} Non-centrality parameter of the chisq test of association between polygenic score and target trait.}
#'	\item{\code{p} Expected P-value of the chisq test of association between polygenic score and target trait.}
#'	\item{\code{power} Power of the chisq test of association between polygenic score and target trait.}
#'	\item{\code{FDR} Expected proportion of false positives among selected markers.}
#'	\item{\code{AUC} For binary traits, area under ROC curve.}
#'	\item{\code{MSE} For quantitative traits, mean square error between target trait and polygenic score.}
#'	\item{\code{NRI} Net reclassification improvement in cases, controls, and combined.}
#'	\item{\code{IDI} Integrated discrimination improvement.}
#'	\item{\code{error} Error message, if any.}
#' }

#' @examples
#' # P-value for ISC schizophrenia score associated with schizophrenia in MGS-EA
#' # See page 3, column 2, paragraph 3 of Dudbridge (2013)
#' polygenescore(74062,n=c(3322+3587,2687+2656),vg1=0.269,pi0=0.99,binary=TRUE,
#' sampling=c(3322/6909,2687/5343),pupper=c(0,0.5),prevalence=.01)$p
#' # [1] 1.029771e-28
#'
#' # Power for ISC schizophrenia score associated with bipolar disorder in WTCCC
#' # See page 4, column 2, paragraph 2 of Dudbridge (2013)
#' polygenescore(74062,c(3322+3587,1829+2935),vg1=0.287,cov12=0.28*0.287,binary=TRUE,
#' sampling=c(3322/6909,1829/4764),pupper=c(0,0.5),prevalence=.01)$power
#' # [1] 0.8042843
#'
#' # Power for cross validation study of Framingham risk score
#' # See page 6, column 1, paragraph 1 of Dudbridge (2013)
#' polygenescore(100000,c(1575,175),vg1=1,pupper=c(0,0.1,0.2,0.3,0.4,0.5),
#' nested=FALSE)$power
#' # [1] 0.19723400 0.11733175 0.09195134 0.07733049 0.06771049
#'
#' # Net reclassification index for cardiovascular disease with QRISK-2 and 53 SNPs
#' # See table 3, row 1, columns 5-6 of Dudbridge et al (submitted)
#' # results vary due to stochastic evaluation of multivariate normal probabilities
#' polygenescore(nsnp=1e5,n=63746+130681,vg1=0.3,pi0=0.8,binary=TRUE,
#' prevalence=0.15,sampling=63746/194427,pupper=c(0,5e-8),
#' r2gx=0.3,r2xy=0.052,corgx=0.1,riskthresh=0.1,adjustedEffects=TRUE)$NRI
#' # [1] -0.006042718  0.015266759  0.009224041

#' @author Frank Dudbridge

#' @references
#' Dudbridge F (2013) Power and predictive accuracy of polygenic risk scores. PLoS Genet 9:e1003348
#' @references
#' Dudbridge F, Pashayan N, Yang J.  Predictive accuracy of combined genetic and environmental risk scores.  Submitted.
#' @export
polygenescore=function(nsnp,
                       n,
                       vg1=0,
                       cov12=vg1,
                       pi0=0,
                       pupper=c(0,1),
		            nested=TRUE,
                       weighted=TRUE,
                       binary=c(FALSE,FALSE),
                       prevalence=c(0.1,0.1),
                       sampling=prevalence,
                       lambdaS=NA,
                       shrinkage=FALSE,
                       logrisk=FALSE,
                       alpha=0.05,
		            r2gx=0,
		            corgx=0,
		            r2xy=0,
  		 		 adjustedEffects=FALSE,
		     		 riskthresh=0.1) {

errorMsg=""
if (min(pupper)<0 | max(pupper)>1)
  errorMsg="Error: entries in pupper should be in (0,1)"
if (errorMsg!="")
  return(list(R2=NULL,NCP=NULL,p=NULL,power=NULL,AUC=NULL,MSE=NULL,error=errorMsg))

if (length(n)==1) n=rep(n,2)
if (length(binary)==1) binary=rep(binary,2)
if (length(prevalence)==1) prevalence=rep(prevalence,2)
if (length(sampling)==1) sampling=rep(sampling,2)

# if including environmental score, force the target trait to be the same as training trait
if (r2xy>0) {
  binary=rep(binary[1],2)
}

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

# covariance of gene score with X
  covgx = corgx*sqrt(vg1*r2gx*r2xy)

# variance of nonzero betas
betaVar = vg1/(nsnp*(1-pi0))
if (r2xy==0) betaCov = cov12/(nsnp*(1-pi0))
else betaCov = betaVar
betaCovX = covgx/(nsnp*(1-pi0))

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
if (binary[1]) betaCovX = betaCovX*liab2obs[1]^2

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
# regression coefficient in SNPs with effects on training trait
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

# covariance between environmental score and estimated gene score
  if (weighted) {
# regression coefficient in SNPs with effects on X
    scoreCovarianceX = betaCovX/(betaVar+samplingVar)
# covariance in SNPs with effects
    scoreCovarianceX = scoreCovarianceX*varBetaHat*probTruncBeta
  }
  else {
    scoreCovarianceX=rep(0,length(plower))
    for(i in 1:length(plower)) {
      scoreCovarianceX[i] = 2*betaCovX/betaVar*(1-pi0)*nsnp*
        integrate(discordantSign,0,Inf,sqrt(betaVar),betaHatThreshLo[i],betaHatThreshHi[i],sqrt(samplingVar),abs.tol=1e-12,stop.on.error=F)$value
    }
  }

# add in effect of environmental risk score
  if (r2xy>0) {
    if (r2gx>1 | r2gx<0) {
      mynan=rep(NaN,length(plower))
      return(list(error="r2gx should be between 0 and 1"))
    }

    cbeta=c(1,1)
    if (adjustedEffects) {
      sigma=matrix(c(1,scoreCovarianceX/r2xy,scoreCovarianceX/varGeneScoreHat,1),nrow=2)
      if (binary[1]) sigma[2,1]=sigma[2,1]/liab2obs[1]^2
      cbeta=solve(sigma)%*%c(scoreCovariance/varGeneScoreHat,1)
    }
    scoreCovariance = cbeta[1]*scoreCovariance+cbeta[2]*r2xy
    if (binary[1] && binary[2]) scoreCovariance = scoreCovariance+cbeta[2]*r2xy*(liab2obs[1]*liab2obs[2]-1)
    if (binary[1] && !binary[2]) scoreCovariance = scoreCovariance+cbeta[2]*r2xy*(liab2obs[1]-1)
    if (!binary[1] && binary[2]) scoreCovariance = scoreCovariance+cbeta[2]*r2xy*(liab2obs[2]-1)
    varGeneScoreHat = cbeta[1]^2*varGeneScoreHat+cbeta[2]^2*r2xy+2*cbeta[1]*cbeta[2]*scoreCovarianceX
    if (binary[1]) varGeneScoreHat = varGeneScoreHat+cbeta[2]^2*r2xy*(liab2obs[1]^2-1)
  }


# Coefficient of determination!
  R2 = scoreCovariance^2/varGeneScoreHat/varY2

# catch impossible configurations
  if (TRUE) {
    if (is.nan(max(R2))) { # } | max(R2)>1 | vg1>0 & abs(cov12)/sqrt(vg1)>1) {
      mynan=rep(NaN,length(plower))
      return(list(R2=mynan,NCP=mynan,p=mynan,power=mynan,AUC=mynan,MSE=mynan,error="cov12 incompatible with vg1"))
    }
  }

# Non-centrality parameter!
  NCP=n[2]*R2/(1-R2)
  if (max(R2)>=1) {
    NCP[R2>=1]=1e10
    errorMsg="Warning: R2>=1 on observed scale"
  }

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
  if (binary[2]) R2 = R2*varY2/liab2obs[2]^2

# net reclassifation improvement!
  if (r2xy>0 & binary[1]) {
# risk thresholds for existing (X) and new (S) scores
    if (logrisk) {
      nriThreshS=log(riskthresh)
      nriThreshX=log(riskthresh)
    }
    else {
      nriThreshS=qnorm(riskthresh,sd=sqrt(1-R2))-qnorm(prevalence[2])
      nriThreshX=qnorm(riskthresh,sd=sqrt(1-r2xy))-qnorm(prevalence[2])
    }

# covariance on risk/liability scale between existing and new scores
    covXS = (cbeta[1]*scoreCovarianceX/liab2obs[2]^2+cbeta[2]*r2xy)

# NRI elements
    if (logrisk) {
      mu = c(log(prevalence[2])-varGeneScoreHat/liab2obs[2]^2/2,log(prevalence[2])-r2xy/2)
      # integrated discrimination improvement
      IDI = exp(mu[1]+varGeneScoreHat/liab2obs[2]^2*3/2)-exp(mu[2]+r2xy*3/2)-exp(mu[1]+varGeneScoreHat/liab2obs[2]^2/2)+exp(mu[2]+r2xy/2)
      if (riskthresh==0) { # continuous NRI
        case10 = 2*pnorm(0,mean=mu[2]+r2xy-mu[1]-varGeneScoreHat/liab2obs[2]^2,sd=sqrt(varGeneScoreHat/liab2obs[2]^2+r2xy-2*covXS))-1
	case01 = 0
	control10 = 2*pnorm(0,mean=mu[1]-mu[2],sd=sqrt(varGeneScoreHat/liab2obs[2]^2+r2xy-2*covXS))-1
	control01 = 0
      }
      else { # threshold NRI
        sigma = rbind(c(varGeneScoreHat/liab2obs[2]^2,covXS),c(covXS,r2xy))
        case10 = mvtnorm::pmvnorm(lower=c(nriThreshS,-Inf),upper=c(Inf,nriThreshX),mean=mu+c(varGeneScoreHat/liab2obs[2]^2,r2xy),sigma=sigma)
        case01 = mvtnorm::pmvnorm(lower=c(-Inf,nriThreshX),upper=c(nriThreshS,Inf),mean=mu+c(varGeneScoreHat/liab2obs[2]^2,r2xy),sigma=sigma)
        control10 = mvtnorm::pmvnorm(lower=c(-Inf,nriThreshX),upper=c(nriThreshS,Inf),mean=mu,sigma=sigma)
        control01 = mvtnorm::pmvnorm(lower=c(nriThreshS,-Inf),upper=c(Inf,nriThreshX),mean=mu,sigma=sigma)
      }
    }
    else {
      diseasethresh=qnorm(prevalence[2],lower=F)
      # integrated discrimination improvement
      thresholdDensity = dnorm(qnorm(prevalence[2]))/prevalence[2]
      oldCaseMean = thresholdDensity*r2xy
      oldCaseVariance = r2xy*(1-oldCaseMean*(thresholdDensity+qnorm(prevalence[2])))
      thresholdDensity = dnorm(qnorm(prevalence[2]))/(1-prevalence[2])
      oldControlMean = -thresholdDensity*r2xy
      oldControlVariance = r2xy*(1+oldControlMean*(thresholdDensity-qnorm(prevalence[2])))
      # approx IDI: risk of the expected liability
      #IDI = pnorm(diseasethresh,caseMean,sd=sqrt(1-caseVariance),lower=F)
      #IDI = IDI - pnorm(diseasethresh,oldCaseMean,sd=sqrt(1-oldCaseVariance),lower=F)
      #IDI = IDI - pnorm(diseasethresh,controlMean,sd=sqrt(1-controlVariance),lower=F)
      #IDI = IDI + pnorm(diseasethresh,oldControlMean,sd=sqrt(1-oldControlVariance),lower=F)
      #
      # exact IDI: expected risk of the liability
      risk=function(l) {
        pnorm(diseasethresh,l,sqrt(1-r2xy))*(dnorm(l,oldCaseMean,sqrt(oldCaseVariance))-
                                             dnorm(l,oldControlMean,sqrt(oldControlVariance)))-
        pnorm(diseasethresh,l,sqrt(1-R2))*(dnorm(l,caseMean,sqrt(caseVariance))-
                                           dnorm(l,controlMean,sqrt(controlVariance)))
      }
      IDI = integrate(risk,-Inf,Inf)$value
      gamma = sqrt(R2/varGeneScoreHat*liab2obs[2]^2)
      if (riskthresh==0) { #continuous NRI
        sigma = rbind(c(1,R2-r2xy),c(R2-r2xy,R2+r2xy-2*gamma*covXS))
        case10 = 2*mvtnorm::pmvnorm(lower=c(qnorm(prevalence[2],lower=F),0),upper=c(Inf,Inf),sigma=sigma)/prevalence[2]-1
     	  case01 = 0
	  control10 = 2*mvtnorm::pmvnorm(lower=c(-Inf,-Inf),upper=c(qnorm(prevalence[2],lower=F),0),sigma=sigma)/(1-prevalence[2])-1
	  control01 = 0
      }
      else { # threshold NRI
        sigma = rbind(c(1,R2,r2xy),c(R2,R2,gamma*covXS),c(r2xy,gamma*covXS,r2xy))
        case10 = mvtnorm::pmvnorm(lower=c(qnorm(prevalence[2],lower=F),nriThreshS,-Inf),upper=c(Inf,Inf,nriThreshX),sigma=sigma)/prevalence[2]
        case01 = mvtnorm::pmvnorm(lower=c(qnorm(prevalence[2],lower=F),-Inf,nriThreshX),upper=c(Inf,nriThreshS,Inf),sigma=sigma)/prevalence[2]
        control10 = mvtnorm::pmvnorm(lower=c(-Inf,-Inf,nriThreshX),upper=c(qnorm(prevalence[2],lower=F),nriThreshS,Inf),sigma=sigma)/(1-prevalence[2])
        control01 = mvtnorm::pmvnorm(lower=c(-Inf,nriThreshS,-Inf),upper=c(qnorm(prevalence[2],lower=F),Inf,nriThreshX),sigma=sigma)/(1-prevalence[2])
      }
    }

    NRIcase = as.numeric(case10-case01)
    NRIcontrol = as.numeric(control10-control01)
  }

  else {
    NRIcase = NULL
    NRIcontrol = NULL
    IDI = NULL
  }

return(list(R2=R2,NCP=NCP,p=pchisq(NCP+1,1,lower=F),power=power,FDR=probTruncNull/(probTruncBeta+probTruncNull),AUC=AUC,MSE=MSE,NRI=c(NRIcase,NRIcontrol,NRIcase+NRIcontrol),IDI=IDI,error=errorMsg))
}

# numerical integration in covariance for unweighted score
discordantSign=function(x,xsigma,threshLo,threshHi,asigma) {
  x*dnorm(x,sd=xsigma)*
  (pnorm(threshLo,mean=x,sd=asigma)-pnorm(threshHi,mean=x,sd=asigma)-pnorm(-threshHi,mean=x,sd=asigma)+pnorm(-threshLo,mean=x,sd=asigma))

}
#########
#end of polygenescore
#########


