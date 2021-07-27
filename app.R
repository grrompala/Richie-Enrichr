## Load necessary packages
library(shiny)
library(dplyr)
library(enrichR)
library(DT)
library(shinydashboard)
library(shinyIncubator)
library(shinyWidgets)
library(ggplot2)

 
options(shiny.maxRequestSize=30*1024^2)  # max-csv upload set to 30 MB
dbs <- listEnrichrDbs() # Object with all Enrichr Database IDs

############ USER INTERFACE ##########

ui <- dashboardPage(skin="red",
                               
                 
#### Sidebar header

        dashboardHeader(title = span("Richie Enrichr: Tool for Interfacing with
                                Enrichr Web Server for Gene Set Enrichment Analysis", 
                        style = "font-size:16px"),titleWidth = 700
        ),

#### Sidebar Features

        dashboardSidebar(
         
          width=300,
  
            
# Input CSV             
            
            fileInput(inputId = "file",
                    label = "Upload modules as .csv",
                    accept = c(".csv"),placeholder = "Submit CSV to start"),
                    
          
# Select Enrichr Database
            
            selectInput('database', label = 'Select Enrichr database(s)', 
                       choices = dbs,multiple=T),
          
# Set Rank in Single DB
          
            sliderInput(inputId = "rank",
                  label="Set Single-DB # of ranked terms",
                  min=1,
                  max=10,
                  value=1),
          
# Set adj.p cut-off in Multi-DB
          
            sliderInput(inputId = "cutoff",
                  label="Adj. p val for Multi-DB function",
                  min=0,
                  step=0.001,
                  max=1.001,
                  value=0.05),
      
# Action Buttons
          
            actionButton("do", "Start Single-DB"),
            actionButton("two","Start Multi-DB"),
            br(),br(),
          # Download Buttons
            div(downloadButton("downloadData", "Download Single-DB Results",class="fines"),style="text-align:center"),
            tags$head(tags$style(".fines{color: black !important;}")),
            br(),
            br(),
            div(downloadButton("downloadData.two","Download Multi-DB Results",class="fine"),style="text-align:center"),
            tags$head(tags$style(".fine{color: black !important;}"))
            
        ),

#### Main Panels    
        dashboardBody(
           
        
           tabsetPanel(
             
              tabPanel("Instructions",tags$img(src="Example.input.png",height=600,width=900)),
             
              tabPanel("Your Modules",DT::dataTableOutput("Modules"),br(),textOutput("Mess")),
             
              tabPanel("Single-DB function",div(DT::dataTableOutput("GO"),
                       style="font-size: 75%; width: 50%")),
             
              tabPanel("Multi-DB function",DT::dataTableOutput("Multi")),
              
              tabPanel("Plots",plotOutput("Bubble"),plotOutput("Ontology"))
            )
      
         )
)

########### SERVER ###############

server <- function(input,output) {

  # Set CSV Input to Data Frame
  
modules <- reactive({
  req(input$file)
  read.csv(input$file$datapath,header=T)
  })

  # Preview Modules after entering CSV

output$Modules <- 
  
  DT::renderDataTable({head(modules(),n=10)})

output$Mess <- renderText({"Upload CSV to see your gene list here"})                      
  #output$Mess <- renderText("ok")
  

  
  # Single DB Function
  
  Out <- eventReactive(input$do,{
   
         withProgress(message = 'Hold your horses', value = 0,{
        
         modules <- read.csv(input$file$datapath,header=T)
       
         all.modules <- unique(modules$Modules)
       
         Output <- data.frame()
         for(x in all.modules){
           
            genes <- modules %>% filter(Modules==x) %>% dplyr::select(id)
            gene.id <- genes$id
            output <- enrichr(gene.id, input$database)
            Temp <- output[[1]] %>% arrange(P.value)
            Temp <- Temp[1:input$rank,]
            Temp$module <- x
            Temp$module.size <- length(gene.id)
            range <- 1:10
            Temp$rank <- as.character(range[1:input$rank])
            Output <- rbind(Output,Temp)
            Output <- Output %>% mutate(SIG=-log10(P.value))
          } 
        })
        Output %>% arrange(SIG)
        }) 
		
   # Multiple DB Function
   
  Multi <- eventReactive(input$two,{
    
           withProgress(message = 'Hold your horses', value = 0,{
           
           modules <- read.csv(input$file$datapath,header=T)
           
           all.modules <- unique(modules$Modules)
           
           summary <- data.frame()
 
           for(x in all.modules){
              genes <- modules %>% filter(Modules==x) %>% dplyr::select(id)
              gene.id <- genes$id 
              top.dbs <- vector()
  
              for(y in input$database){
                 output <- enrichr(gene.id,y)
                 Temp <- output[[1]] %>% arrange(Adjusted.P.value)
                 Temp.term <- Temp$Term[1]
                 Temp.term <- ifelse(Temp$Adjusted.P.value[1]<input$cutoff,Temp.term,"Non.Sig")
                 top.dbs <- append(top.dbs,Temp.term)
                 }
  
           mod.summary <- append(x,top.dbs)
           summary <- rbind(summary,mod.summary)
           colnames(summary) <- c("Module.ID",input$database)
           }
           }) 
           summary
           })

  
 
  # Reactive output for single db function
  observeEvent(input$do,{ 
     
    options(digits = 3)
    output$GO <- DT::renderDataTable({datatable(Out()[,c(10:12,1:4,7:8)])%>%
                 formatRound(c(6:9),3)})
    
  })
  
  
 # Bubble plot with ggplot2
  
  output$Bubble <- renderPlot({
    
    ggplot(Out(),aes(x=Odds.Ratio,-log10(P.value),size=module.size,label=Term))+
     geom_point(alpha=.25,shape=21,colour="black",fill="blue",stroke=2)+
     geom_text(check_overlap=FALSE,aes(label=ifelse(-log10(P.value)>3,Term,""),size=1))+ #,fontface="bold"))+ #hjust=ifelse(pval>3,.1,.3315),vjust= ifelse(OR<60,.6,-.4),fontface="bold"))+
      scale_size_continuous(range = c(3, 20),breaks=c(10,100,1000))+
    # xlim(0,20)+
      ylab("-log10(Adjusted.P.value)")+
      xlab("Odds Ratio")+
     labs(size="Module Size")+
      theme(
       panel.background = element_rect(colour="black",fill="white"),
       legend.title=element_text(size=20),
        legend.text=element_text(size=18),
        axis.text.x=element_text(size=22,face="bold"),
        axis.text.y=element_text(size=22,face="bold"),
        axis.title.x=element_text(size=22,face="bold"),
        axis.title.y = element_text(size=22,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)
      )
})
  
 
  
  output$Ontology <- renderPlot({
    ggplot(Out(),aes(x=SIG,y=reorder(Term,SIG)))+
      geom_bar(stat="identity",fill="blue")+
      xlab("-log10(P.value)")+
      ylab("Term")+
      theme(
        axis.title.x=element_text(size=22,face="bold"),
        axis.text.y=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=22,face="bold"),
        axis.text.x=element_text(size=16,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)
      )
    
  })
      #  scale_fill_manual(name="Mileage", 
      #                   labels = c("Above Average", "Below Average"), 
      #                  values = c("above"="#00ba38", "below"="#f8766d")) + 
      #labs(subtitle="Normalised mileage from 'mtcars'", 
      #       title= "Diverging Bars") + 
      #  coord_flip()
  
  
  
  
   
  # Reactive output for multi-db function
  observeEvent(input$two,{
    
    output$Multi <- DT::renderDataTable({datatable(Multi())}) 
  
  })
  
# Downloadable csvs of selected dataset ----
  
  # Single DB output
  output$downloadData <- 
  
    downloadHandler(
      filename = function() {
      paste("report", ".csv", sep = "")
    },
      content = function(file) {
      write.csv(Out(), file, row.names = TRUE)
    }
    )
  
  # Multi DB output
	output$downloadData.two <- 
	  
	  downloadHandler(
      filename = function() {
      paste("report", ".csv", sep = "")
    },
      content = function(file) {
      write.csv(Multi(), file, row.names = TRUE)
      }
    )

  
  
 }
  

shinyApp(ui, server)