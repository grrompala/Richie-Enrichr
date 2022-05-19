## Load necessary packages
library(shiny)
library(dplyr)
library(enrichR)
library(DT)
library(shinydashboard)
library(shinyIncubator)
library(shinyWidgets)
library(ggplot2)
library(plotly)
library(shinyjs)
library(shinycssloaders)

 
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
                       choices = dbs$libraryName,selected=NULL,multiple=T),
          
# Set Rank in Single DB
 
              sliderInput(inputId = "rank",
                  label="Max number of significant terms to return (per database)",
                  min=1,
                  max=10,
                  value=1,width="100%"),
             
br(),


          
# Set adj.p cut-off in Multi-DB
          
            # sliderInput(inputId = "cutoff",
            #       label="Adj. p val for Multi-DB function",
            #       min=0,
            #       step=0.001,
            #       max=1.001,
            #       value=0.05),
      
# Action Buttons
          
            div(style="padding-left: 20px;",actionBttn("run","Run Enrichr Query",color = "danger",style = "bordered",size = "lg")),
           # tags$head(tags$style(".run{text-align: center !important;}")),
            br(),br(),
          # Download Buttons
          #  div(downloadButton("downloadData", "Download Single-DB Results",class="fines"),style="text-align:center"),
          #  tags$head(tags$style(".fines{color: black !important;}")),
          #  br(),
           # br(),
            div(downloadButton("downloadData.two","Download Results",class="fine"),style="text-align:center"),
            tags$head(tags$style(".fine{color: black !important;}"))
),
            

#### Main Panels    
        dashboardBody(
           
          useShinyjs(),
        
           tabsetPanel(id = "Tabs",
             
              tabPanel("Introduction",
                       
                       
                       h2("Welcome to the Richie Enrichr App"),br(),
                       
                       HTML("<h4><b>Purpose</b>: Perform exploratory enrichment analyses on multiple gene lists in Enrichr Web App<br><br>
                            <ul><li><a href='https://maayanlab.cloud/Enrichr/' target='_blank'>Enrichr</a> is a web app developed
                            by the Ma'ayan lab in the Department of Computational Genomics at Icahn School of Medicine at Mount Sinai.</li><br>
                            <li>The Enrichr app does not take multiple gene lists as input nor make it convenient to download outputs across
                            multiple databases: this app provides that functionality</li></ul></h4>"),
                       br(),
                       fluidRow(
                         box(width = 4,solidHeader = TRUE,status="danger",
                           HTML("Screenshot of <a href='https://maayanlab.cloud/Enrichr/' target='_blank'>Enrichr</a> App"),
                           br(),br(),
                           tags$img(src="Enrichr.png",width="100%")
                           ),
                         box(width=8,title="Instructions",solidHeader = TRUE,status="danger",
                           tags$img(src="Example.input.png",width="100%")
                         )
                       )
                       ),
             
              tabPanel("Your Modules",
                       
                       
                       plotOutput("modsummary_plot") %>% withSpinner(color = "#FF0000"),
                       
                       br(),
                       
                       DT::dataTableOutput("Modules") %>% withSpinner(color = "#FF0000"),
                       
                       br(),
  
                       textOutput("Mess")),
             
             
              tabPanel("Table Outputs",DT::dataTableOutput("Multi") %>% withSpinner(color = "#FF0000")),
              
              tabPanel("Plots",
                       fluidPage(
                       br(),
                       
                       uiOutput("bubble_databases"),
                       
                      # textOutput("clickevent"),
                       
                       box(title="Top Terms for each Module",width=12,status="danger",solidHeader = TRUE,
                       plotlyOutput("Bubble") %>% withSpinner(color = "#FF0000"),
                       DTOutput("bubble_table")
                       ),
                           
                      # DTOutput("bubbleDT"),
                      
                       
                       
                       plotOutput("Ontology")
                      )
              )
            )
      
         )
)

########### SERVER ###############

server <- function(input,output,session) {
  
  
  
  # what tabs to show
  
  observe({
    if(is.null(input$file)){
      hideTab("Tabs","Your Modules")
      
    }else{showTab("Tabs","Your Modules")
    }
    
    if(input$run<1){
      hideTab("Tabs","Table Outputs")
    }else{showTab("Tabs","Table Outputs")}
    
    if(input$run<1){
      hideTab("Tabs","Plots")
    }else{showTab("Tabs","Plots")}
    
    if(is.null(input$file) | is.null(input$database)){
      shinyjs::disable("run")
    }else{
      shinyjs::enable("run")
    }

    
    
  })
  

  # Set CSV Input to Data Frame
  
modules <- reactive({
  req(input$file)
  read.csv(input$file$datapath,header=T)
  })

  # Preview Modules after entering CSV

output$Modules <- 
  
  DT::renderDT({
    modules()
    },filter="top",rownames=FALSE,options=list(pageLength=100)
  )

output$Mess <- renderText({
  if(is.null(modules())==F){
  "Upload CSV to see your gene list here"
    }
  })                      
  #output$Mess <- renderText("ok")
  

  # Module summary
output$modsummary_plot <-
 renderPlot({
  sum <- modules() %>%
     count(Modules)
   
   ggplot(sum, aes(n))+
     geom_histogram(binwidth = 1)+
     xlab("Module Size")+
     ylab("Number of Modules")

 })


  
  # Single DB Function
  
  Out <- eventReactive(input$do,{
   
         withProgress(message = 'Please wait', value = 0,{
        
         modules <- read.csv(input$file$datapath,header=T)
       
         all.modules <- unique(modules$Modules)
       
         Output <- data.frame()
         for(x in all.modules){
           
            genes <- modules %>% filter(Modules==x) %>% dplyr::select(id)
            gene.id <- genes$id
            output <- enrichr(gene.id, input$database)[[1]]
            if(nrow(output)>0){
            Temp <- output %>% arrange(P.value)
            Temp <- Temp[1:input$rank,]
            Temp$module <- x
            Temp$module.size <- length(gene.id)
            range <- 1:10
            Temp$rank <- as.character(range[1:input$rank])
            Output <- rbind(Output,Temp)
          } 
         }
         })
        Output %>% mutate(SIG=-log10(P.value)) %>% arrange(SIG)
        }) 
		
   # Multiple DB Function
   
  Multi <- eventReactive(input$run,{
    
           withProgress(message = 'Please wait', value = 4,{
           
           modules <- read.csv(input$file$datapath,header=T)
           
           all.modules <- unique(modules$Modules)
           
 
           lapply(all.modules,function(x){
              genes <- modules %>% filter(Modules==x) %>% dplyr::pull(id)
              
           lapply(input$database,function(db){
                 output <- enrichr(genes,db)
                 Temp <- output[[1]] %>% arrange(Adjusted.P.value) %>% .[1:input$rank,]
                 Temp %>% mutate(Database=db)
                 }) %>% bind_rows() %>% mutate(Module_ID=x) %>% mutate(Module.size=length(genes))

           }) %>% bind_rows()
           
    
           })
  })

  
 
  # Reactive output for single db function
  # observeEvent(input$do,{ 
  # 
  #   output$GO <- DT::renderDataTable({
  #     
  #     datatable(Out() %>% select(module,
  #                                module.size,
  #                                rank,
  #                                "Name"=Term,
  #                                Odds.Ratio,
  #                                P.value,
  #                                Adjusted.P.value,
  #                                Overlap,
  #                                Genes)
  #     ) %>% formatRound(c(5:7),3)
  #       })
  #   
  # })
  
  
 # Bubble plot with ggplot2
  
  output$bubble_databases <- renderUI({
    req(is.null(input$database)==F)
    selectInput("bubble_db","Choose Database or use top term across all",choices=c("all",input$database),selected = "all")
  })
  
#  output$clickevent <- renderPrint({
#   event_data("plotly_click")
#  })
  
  
  bubble_data <- reactive({Multi() %>% {if(input$bubble_db !="all") dplyr::filter(.,Database==input$bubble_db) else .} %>% group_by(Module_ID) %>%
    dplyr::filter(Adjusted.P.value==min(Adjusted.P.value)) %>% distinct(Module_ID,Adjusted.P.value,.keep_all=TRUE)
  })
  
  output$Bubble <- renderPlotly({
    

    
    fig <- plot_ly(bubble_data(), x = ~Odds.Ratio,
                   y = ~-log10(Adjusted.P.value),
                   text = ~paste("Module ID:",Module_ID,"\n",Database,":",Term),
                   size = ~Module.size,
                   type = 'scatter',
                   mode = 'markers',
                   marker = list(symbol = 'circle', sizemode = 'diameter',
                                 opacity = 0.5,
                                 line = list(width = 2, color = '#000000')),
                   sizes = c(10, 50),
                   color=~Database)
    fig <- fig %>% layout(
                          xaxis = list(showgrid = FALSE),
                          yaxis = list(showgrid = FALSE))
    
    fig
    
  })
    
  observeEvent(event_data("plotly_click"),{
    output$bubble_table <- renderDT({
   
    cords <- event_data("plotly_click")
    
    curve <- sort(input$database)[cords$curveNumber+1] 
    
      
    mod <- bubble_data() %>% {if(input$bubble_db=="all") dplyr::filter(.,Database==curve) else .} %>% .[cords$pointNumber+1,] %>% pull(Module_ID)
  

    #  mod <- Multi() %>%
    #    filter(round(Odds.Ratio,digits = 1)==cords$x) %>%
    #    pull(Module_ID)
      
      Multi() %>% filter(Module_ID==mod)
    },rownames=FALSE,options=list(scrollX = TRUE))
  })
    
    
    
    # ggplot(data,aes(x=Odds.Ratio,-log10(Adjusted.P.value),size=Module.size,label=Term))+
    #  geom_point(alpha=.25,shape=21,colour="black",fill="blue",stroke=2)+
    #  geom_text(check_overlap=FALSE,aes(label=ifelse(-log10(Adjusted.P.value)>1.3,Term,""),size=1))+ #,fontface="bold"))+ #hjust=ifelse(pval>3,.1,.3315),vjust= ifelse(OR<60,.6,-.4),fontface="bold"))+
    #   scale_size_continuous(range = c(3, 20),breaks=c(10,100,1000))+
    # # xlim(0,20)+
    #   ylab("-log10(Adjusted.P.value)")+
    #   xlab("Odds Ratio")+
    #  labs(size="Module Size")+
    #   theme(
    #    panel.background = element_rect(colour="black",fill="white"),
    #    legend.title=element_text(size=20),
    #     legend.text=element_text(size=18),
    #     axis.text.x=element_text(size=22,face="bold"),
    #     axis.text.y=element_text(size=22,face="bold"),
    #     axis.title.x=element_text(size=22,face="bold"),
    #     axis.title.y = element_text(size=22,face="bold"),
    #     panel.border = element_rect(colour = "black", fill=NA, size=2)
    #   )

#  output$bubbleDT <- renderDT({bubble_data()},options=list(scrollX = TRUE))
 
  
  output$Ontology <- renderPlotly({
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
  observeEvent(input$run,{
    
    output$Multi <- DT::renderDT({datatable(Multi(),rownames=FALSE)}) 
    
    updateTabsetPanel(session, inputId="Tabs", selected="Plots")
  
  })
  
  observeEvent(input$file,ignoreInit = TRUE,{
  updateTabsetPanel(session, inputId="Tabs", selected="Your Modules")
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