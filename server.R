library(shiny)

shinyServer(function(input, output) {
  source("avengeme.R")

  results = eventReactive(input$go, {

    p=as.numeric(strsplit(input$p,",")[[1]])
    pupper=as.numeric(strsplit(input$pupper,",")[[1]])
    vg1=input$vg1
    if (input$estimateVariance) vg1=NA
    vg2=input$vg2
    if (input$bidirectional & input$estimateVg2) vg2=NA
    cov12=input$cov12
    if (input$estimateCovariance) cov12=NA
    pi0=input$pi0
    if (input$estimatePi0) pi0=NA
    pi02=input$pi02
    if (input$bidirectional & input$estimatePi02) pi02=NA
    if (input$fix) {
      vg2=vg1
      cov12=vg1
      pi02=pi0
    }

    if (!is.na(vg1) & !is.na(cov12) & !is.na(pi0) & (!input$bidirectional | !is.na(vg2) & !is.na(pi02))) {
      result=polygenescore(nsnp=input$nsnp,n=c(input$n1,input$n2),vg=vg1,cov12=cov12,pi0=pi0,pupper=pupper,nested=input$nested,weighted=input$weighted,binary=c(input$binary1,input$binary2),prevalence=c(input$prevalence1,input$prevalence2),sampling=c(input$sampling1,input$sampling2))

      outputString="<HR><H4>No parameters estimated</H4>"
      if (result$error=="") {
        outputString=paste(outputString,"R2 between score and target trait = ",toString(format(result$R2,digits=3)))
        outputString=paste(outputString,"<BR>Non-centrality parameter of chisq test = ",toString(format(result$NCP,digits=3)))
        outputString=paste(outputString,"<BR>Expected p-value of chisq test = ",toString(format(result$p,digits=3)))
        outputString=paste(outputString,"<BR>Power of chisq test = ",toString(format(result$power,digits=3)))
        if (input$binary2) {
          outputString=paste(outputString,"<BR>Area under ROC curve = ",toString(format(result$AUC,digits=3)))
        }
        else {
          outputString=paste(outputString,"<BR>Mean square error = ",toString(format(result$MSE,digits=3)))
        }
      }
      else {
        outputString=paste(outputString,"<font color=#ff0000>",result$error,"</font>")        
      }
      outputString
    }
    else {
      result=estimatePolygenicModel(p=p,nsnp=input$nsnp,n=c(input$n1,input$n2),vg=c(vg1,vg2),cov12=cov12,pi0=c(pi0,pi02),pupper=pupper,nested=input$nested,weighted=input$weighted,binary=c(input$binary1,input$binary2),prevalence=c(input$prevalence1,input$prevalence2),sampling=c(input$sampling1,input$sampling2),bidirectional=input$bidirectional,fixvg2pi02=input$fix)

      outputString="<HR>"
      if (result$error=="") {
        outputString=paste(outputString,"<H4>Estimated polygenic model</H4>")

        if (input$bidirectional) {
          if (is.na(vg1)) outputString=paste(outputString,"Variance in training sample = ",format(result$vg[1,1],digits=3),"(",format(result$vg[1,2],digits=3),",",format(result$vg[1,3],digits=3),")")
          if (!input$fix & is.na(vg2)) outputString=paste(outputString,"<BR>Variance in target sample = ",format(result$vg[2,1],digits=3),"(",format(result$vg[2,2],digits=3),",",format(result$vg[2,3],digits=3),")")
          if (!input$fix & is.na(cov12)) outputString=paste(outputString,"<BR>Covariance between training and target samples = ",format(result$cov12[1],digits=3),"(",format(result$cov12[2],digits=3),",",format(result$cov12[3],digits=3),")")
          if (is.na(pi0)) outputString=paste(outputString,"<BR>Proportion of null SNPs in training sample = ",format(result$pi0[1,1],digits=3),"(",format(result$pi0[1,2],digits=3),",",format(result$pi0[1,3],digits=3),")")
          if (!input$fix & is.na(pi02)) outputString=paste(outputString,"<BR>Proportion of null SNPs in target sample = ",format(result$pi0[2,1],digits=3),"(",format(result$pi0[2,2],digits=3),",",format(result$pi0[2,3],digits=3),")")
        }

        else {
          if (is.na(vg1)) outputString=paste(outputString,"Variance in training sample = ",format(result$vg[1],digits=3),"(",format(result$vg[2],digits=3),",",format(result$vg[3],digits=3),")")
          if (!input$fix & is.na(cov12)) outputString=paste(outputString,"<BR>Covariance between training and target samples = ",format(result$cov12[1],digits=3),"(",format(result$cov12[2],digits=3),",",format(result$cov12[3],digits=3),")")
          if (is.na(pi0)) outputString=paste(outputString,"<BR>Proportion of null SNPs in training sample = ",format(result$pi0[1],digits=3),"(",format(result$pi0[2],digits=3),",",format(result$pi0[3],digits=3),")")
        }
      outputString=paste(outputString,"<BR>Log-likelihood = ",format(result$logLikelihood,digits=3))
      }

      else {
        outputString=paste(outputString,"<font color=#ff0000>",result$error,"</font>")
      }
      outputString=paste(outputString,"<HR>")
      outputString
    }
  })

  output$polygenescore=renderText({results()})

}
)
