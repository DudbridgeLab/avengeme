library(shiny)

shinyUI(fluidPage(

  titlePanel("AVENGEME"),
  p("Additive Variance Explained and Number of Genetic Effects Method of Estimation"),
  p("Luigi Palla and Frank Dudbridge, Am J Hum Genet (2015) in press"),

  fluidRow(
    column(2,
      numericInput("nsnp",label="Number of SNPs",value=100000),
      numericInput("n1",label="Training sample size",value=1000),
      numericInput("n2",label="Target sample size",value=1000)
),

column(4,
      checkboxInput("fix",label="Identical models in training and target samples",value=FALSE),

      checkboxInput("bidirectional",label="Bidirectional estimation",value=FALSE),

      checkboxInput("estimateVariance","Estimate variance in training sample",value=TRUE),
      conditionalPanel(
        condition = "input.estimateVariance == false",
        sliderInput("vg1",label="Fix variance to value:",value=0.5,min=0,max=1,step=0.01)
      ),
        
      conditionalPanel(
        condition = "input.fix == false",
        checkboxInput("estimateCovariance","Estimate covariance between training and target samples",value=TRUE),
        conditionalPanel(
          condition = "input.estimateCovariance == false",
          sliderInput("cov12",label="Fix covariance to value:",value=0.5,min=0,max=1,step=0.01)
        )
      ),

      conditionalPanel(
        condition = "input.fix == false & input.bidirectional == true",
        checkboxInput("estimateVg2","Estimate variance in target sample",value=TRUE),
        conditionalPanel(
          condition = "input.estimateVg2 == false",
          sliderInput("vg2",label="Fix variance to value:",value=0.5,min=0,max=1,step=0.01)
        )
       ),

      checkboxInput("estimatePi0","Estimate proportion of null SNPs in training sample",value=TRUE),
      conditionalPanel(
        condition = "input.estimatePi0 == false",
        sliderInput("pi0",label="Fix null proportion to value:",value=0.95,min=0,max=1,step=0.01)
      ),

      conditionalPanel(
        condition = "input.fix == false & input.bidirectional == true",
        checkboxInput("estimatePi02","Estimate proportion of null SNPs in target sample",value=TRUE),
        conditionalPanel(
          condition = "input.estimatePi02 == false",
          sliderInput("pi02",label="Fix null proportion to value:",value=0.95,min=0,max=1,step=0.01)
        )
      )
),

column(3,
      checkboxInput("binary1",label="Binary training trait",value=FALSE),
      conditionalPanel(
        condition = "input.binary1 == true",
        sliderInput("prevalence1",label="Training trait prevalence",value=0.5,min=0,max=1,step=0.01),
        sliderInput("sampling1",label="Training sampling fraction",value=0.5,min=0,max=1,step=0.01)
      ),

      checkboxInput("binary2",label="Binary target trait",value=FALSE),
      conditionalPanel(
        condition = "input.binary2 == true",
        sliderInput("prevalence2",label="Target trait prevalence",value=0.5,min=0,max=1,step=0.01),
        sliderInput("sampling2",label="Target sampling fraction",value=0.5,min=0,max=1,step=0.01)
      )

),

column(3,
      textInput("pupper",label="P-value selection thresholds in training sample",value="0,1"),
      checkboxInput("nested",label="Nested intervals",value=TRUE),
      checkboxInput("weighted",label="Weighted score",value=TRUE),

      textInput("p",label="P-values or Z-scores between polygenic score and target trait",value=""),
      actionButton("go","Go! Estimate polygenic model")

)


    ),

    htmlOutput("polygenescore")
)
)


