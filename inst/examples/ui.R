

library(shinythemes)
library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(shinytoastr)

shinyUI(fluidPage(
  theme=shinytheme("yeti"),
  shinyjs::useShinyjs(),
  useToastr(),
  # themeSelector(),
  pageWithSidebar(
    headerPanel(title=HTML("Genomic prediction"),
                windowTitle="Genomic prediction"),
    sidebarPanel(


    ## conditionalPanel() functions for selected tab
    conditionalPanel(condition="input.tabselected==1",
                     h3("Genotypic information"),

                     # h6('follow the step by step'),

                     fileInput(inputId='file1',
                               label="Load marker information in '.in' format",
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv',"in"
                               )),

                     h3("Sample's names "),
                     fileInput(inputId='file2',
                               label="Load file with samples in '.txt' format",
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv'
                               )),

                     h3('Phenotypic information'),
                     fileInput(inputId='file3',
                               label="Load file with phenotypic information in '.csv' format",
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv'
                               )),
                     # checkboxGroupInput("checkGroup", label = h3("Models"),
                     #                    choices = list("ASReml" = 1, "RKHS" = 2, "sommer" = 3,
                     #                                   "BRR" = 4,  "BayesA"= 5, "BayesB"=6,
                     #                                   "BayesC"=7, "BLasso"=8),
                     #                    selected = c(1,2,3),inline = T),
                     awesomeCheckboxGroup(
                       inputId = "checkGroup",
                       label = h3("Models"),
                       choices = c("ASReml","RKHS","sommer","BRR","BayesA","BayesB","BayesC","BLasso"),
                       selected = c("RKHS"),
                       inline = TRUE,
                       status = "danger"
                     ),
                     pickerInput(
                       inputId = "Id094",
                       label = h3("Select the traits"),
                       choices = NULL,
                       options = list(size = 5,#style = "btn-danger",
                         `actions-box` = TRUE),
                       multiple = TRUE
                     ),

                     sliderInput("porcent", label = h4("Percentage of test population"), min = 0,
                                 max = 0.98, value = 0.3, step = 0.1),
                     div(id = "nonIter",
                     sliderInput("iter", label = h4("Number of iterations"), min = 1,
                                 max = 100, value = 2, step = 1, post = " iterations")),
                     actionBttn(
                       inputId = "action",
                       label = "Run !",
                       style = "jelly",
                       color = "warning",icon = icon("sync")
                     )
                     # actionButton( inputId= "action", label = "Run", icon("sync"))
                     ),


    conditionalPanel(condition="input.tabselected==2",
                     radioButtons("choice","Choose an option", choices=c("Dataset" = 1 )),
                     downloadButton("downloadData", "Download results")

    ),


    conditionalPanel(condition="input.tabselected==3",
                     h3("Prediction ability"),
                     ),

    conditionalPanel(condition="input.tabselected==4",
                     h3("Marker information"),
                     awesomeRadio(
                       inputId = "method",
                       label = "Radio buttons",
                       choices = c("BRR","BayesA","BayesB","BayesC","BLasso"),
                       selected = "BRR",
                       inline = TRUE,
                       checkbox = T
                     ),
                     pickerInput(
                       inputId = "Id095",
                       label = "Select the traits",
                       choices = NULL
                     )
    )
    # img(src="CIAT1.png", height = 55, width = 120),
    # img(src="univ.jpg", height = 60, width = 50),
    # br(),
    # tags$em("johan.aparicio@correounivalle.edu.co")

  ),


  mainPanel(

    tabsetPanel(
      tabPanel("Visualization", value=1,
               h5("First load the files, then select run and wait... "),br(),
               fluidRow(column(withSpinner(DT::dataTableOutput("Rawdata"),type = 5,color = "#337ab7"),width = 12) ) ),

      tabPanel("Progress", value=2, conditionalPanel(condition="input.choice==1", br() ,
                                                     div(style = 'text-align: center',pre(id = "console")),
                                                 br() ),
               conditionalPanel(condition="input.choice==2", verbatimTextOutput("struct"))),

      tabPanel("Plot", value=3,br(), plotOutput("plot"),
               br(),plotlyOutput("boxplot"), br()),
      id = "tabselected",
      tabPanel("Marker information",value=4,
               br(),
               h5("If you need consult about how programming this in R, the next link would be helpful"),
               fluidRow(column(withSpinner(plotOutput("mark"),type = 5,color = "#337ab7"),
                               withSpinner(plotOutput("pred"),type = 5,color = "#337ab7"),width = 12) ),
               br()
               )
    )
  )
  ))
)
