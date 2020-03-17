library(shiny)
library(ggplot2)
library(plotly)
library(SpATS)
library(ggsci)
library(ggpubr)
library(shinyjs)
library(shinytoastr)
source("https://raw.githubusercontent.com/AparicioJohan/GPrediction/master/crossGP.R")

options(shiny.maxRequestSize = 70*1024^2)

shinyServer(function(input,output,session)({

# Update variable traits --------------------------------------------------

  observe({
    phen <-  input$file3$datapath
    validate(need(input$file3, "missing file 3"))
    phen <- read.csv(phen)
    updatePickerInput(session, "Id094", choices=names(phen[,-1]))
  })
  
  observe({
    updatePickerInput(session, "Id095", choices=input$Id094)
  })
  
# Read data and run the models --------------------------------------------
  
  df <- eventReactive(input$action,{
    geno <-  input$file1$datapath
    samp <-  input$file2$datapath
    phen <-  input$file3$datapath
    
    
    if(is.null(input$file1)|is.null(input$file2)|is.null(input$file2)) 
      toastr_error( "Missing files..." ,position = "bottom-right")
    
    validate(need(input$file1, "missing file 1"))
    validate(need(input$file2, "missing file 2"))
    validate(need(input$file3, "missing file 3"))
    
    prior <- input$checkGroup
    # prior[prior=="1"] <- "ASReml"
    
    validate(need(input$checkGroup, "Select one model"))
    validate(need(input$Id094, "Select the traits"))
    
    withCallingHandlers({
          results <- crossGP(geno,samp,phen,prior,niter = input$iter, testporc = input$porcent, traits = input$Id094)
      },
      message = function(m) {
        shinyjs::html("console", m$message, TRUE)
      }
    )
      results
    
  })

# Output data -------------------------------------------------------------
  
  output$Rawdata <- DT::renderDataTable({
    req(nrow(df()$data)>=1)
    df()$data %>% 
      DT::datatable(
        extensions = 'Buttons', filter = 'top', selection="multiple",
        options = list(dom = 'lfrtipB', scrollX = TRUE, pageLength = 10,
                       buttons = c('excel', "csv")))
    })
  
# Plot correlations -------------------------------------------------------

  output$plot <- renderPlot({
    if (input$action==0) {return()}
    else {
      validate(need(input$file1, "missing file 1"))
      validate(need(input$file2, "missing file 2"))
      validate(need(input$file3, "missing file 3"))
      if(length(input$checkGroup)>=1){
        g1 <- df()$data %>% 
          ggplot(aes(x=trait,y=corr,fill=prior))+geom_boxplot()+theme_bw(base_size = 13)+
          theme(axis.text.x = element_text(angle = 70,hjust = 1),legend.title = element_blank())+
          ylab("Prediction Ability")+xlab("")+ scale_fill_simpsons()
        
        if(input$iter==1||input$porcent==0){
          g1 <- df()$data %>% 
            ggplot(aes(x=trait,y=corr,color=prior))+geom_point()+theme_bw(base_size = 13)+
            theme(axis.text.x = element_text(angle = 70,hjust = 1),legend.title = element_blank())+
            ylab("Prediction Ability")+xlab("")+ scale_fill_simpsons()
        }
          
        
      } else {
        g1 <- df()$data %>% 
          ggplot(aes(x=trait,y=corr,fill=trait))+geom_boxplot()+theme_bw(base_size = 13)+
          theme(axis.text.x = element_text(angle = 70,hjust = 1),legend.title = element_blank())+
          ylab("Prediction Ability")+xlab("")+ scale_fill_simpsons()
        
        if(input$iter==1||input$porcent==0){
          g1 <- df()$data %>% 
            ggplot(aes(x=trait,y=corr,color=trait))+geom_point()+theme_bw(base_size = 13)+
            theme(axis.text.x = element_text(angle = 70,hjust = 1),legend.title = element_blank())+
            ylab("Prediction Ability")+xlab("")+ scale_fill_simpsons()
        }
        
      }
      isolate(g1)
    }
  })

# Marker information ------------------------------------------------------

  output$mark <- renderPlot({
    if (input$action==0) {return()}
    else {
      validate(need(input$file1, "missing file 1"))
      validate(need(input$file2, "missing file 2"))
      validate(need(input$file3, "missing file 3"))

      req(input$Id095)
      validate(need(input$porcent==0, "Percentage of test population should be zero"))
      validate(need(input$method%in%input$checkGroup, paste(input$method, "is not include in the fitted models")))
      
      out_table <- df()
      
      bHat <- out_table$models[[input$Id095,input$method]]$ETA[[1]]$b
      MarK <- data.frame(Marker=1:length(bHat),bHat = bHat^2)
      g1 <- MarK %>% 
        ggplot(aes(x=Marker, y=bHat)) + 
        geom_point(size=1.5,color="black") + 
        geom_segment(aes(x=Marker, 
                         xend=Marker, 
                         y=0, 
                         yend=bHat), color="red") + 
        labs(title="Marker Effects", 
             subtitle=paste(input$method, "model",sep = " "), 
             caption="generated by: Mr.Gen", y= 'Estimated Squared-Marker Effect') + 
        theme(axis.text.x = element_text(angle=65, vjust=0.6))+theme_bw(base_size = 15)
      
      isolate(g1)
    }
  })
  
  
  output$pred <- renderPlot({
    if (input$action==0) {return()}
    else {
      validate(need(input$file1, "missing file 1"))
      validate(need(input$file2, "missing file 2"))
      validate(need(input$file3, "missing file 3"))
      
      req(input$Id095)
      validate(need(input$porcent==0, " "))
      req(input$method%in%input$checkGroup)
      
      out_table <- df()
      
      bHat <- out_table$models[[input$Id095,input$method]]$ETA[[1]]$b
      MarK <- data.frame(Marker=1:length(bHat),bHat = bHat^2)
      
      Pred <- data.frame( yHat=out_table$models[[input$Id095,input$method]]$yHat,
                          y=out_table$models[[input$Id095,input$method]]$y)  
      g1 <- 
        Pred %>% 
        ggplot(aes(x=yHat,y=y))  +geom_point(alpha=0.2,size=3)+ 
        geom_abline(slope = 1,intercept = 0,linetype=2, color="red",size=1)+
        geom_smooth(method = "lm",formula = y~x, se = F)+
        labs(y="Phenotype",x='Predicted genomic value', title="Predicted genomic values vs phenotypes")+
        stat_cor()+theme_bw(base_size = 15)
      
      isolate(g1)
    }
  })
  

# update number iterations ------------------------------------------------


  observe({
    if (input$porcent != 0) {
      shinyjs::show("nonIter",animType = "fade",anim = TRUE)
      # enable("nonIter")
    } else {
      shinyjs::hide("nonIter",animType = "fade",anim = TRUE)
      # disable("nonIter")
      toastr_info("Percentage of test population equal to zero: Only one iteration selected",position = "bottom-right")
    }
  })
  
  
  observeEvent(input$run,{

    toastr_info(paste(input$Id094,collapse = " "))
      
    })
  
  
}))


