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
library(tidyr)
library(pheatmap)
library(RColorBrewer)
 
options(shiny.maxRequestSize=30*1024^2)  # max-csv upload set to 30 MB
dbs <- listEnrichrDbs() # Object with all Enrichr Database IDs

COLORS <- rownames(brewer.pal.info)

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
              # 
              # sliderInput(inputId = "rank",
              #     label="Max number of significant terms to return (per database)",
              #     min=1,
              #     max=10,
              #     value=1,width="100%"),
             
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
hidden(
            div(id="DOWNLOAD_BUTTON",downloadButton("downloadData.two","Download Results",class="fine"),style="text-align:center"),
            tags$head(tags$style(".fine{color: black !important;}"))
)
),
            

#### Main Panels    
        dashboardBody(
           
          useShinyjs(),
        
           tabsetPanel(id = "Tabs",
             
              tabPanel("Introduction",
                       
                       
                       h2("Welcome to the Richie Enrichr App"),br(),
                       
                       HTML("<h4><b>Purpose</b>: Pull exploratory enrichment analyses from the Enrichr Web App for each unique gene list (module) provided<br><br>
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
                       fluidPage(
                         
                       br(),
                       
                    box(width=12,
                        status = "danger",
                        solidHeader = TRUE,
                        valueBoxOutput("total_modules"),
                        infoBoxOutput("total_genes"),
                        
                        ),

                       fluidRow(column(
                         6,
                         box(
                           width = 12,
                           status = "danger",
                           solidHeader = TRUE,
                           
                           plotOutput("modsummary_plot") %>% withSpinner(color = "#FF0000"),
                           dropdownButton(circle = FALSE,label = "Exclude Modules",icon=icon("network-wired"),status="danger",up = TRUE,
                                    multiInput(
                                      inputId = "select_mods", label = "Select any modules to exclude from analysis",
                                      choices = c("temp"),
                                      options = list(
                                        enable_search = TRUE,
                                        non_selected_header = "All Modules",
                                        selected_header = "Excluded Modules:"
                                      ), width = "350px"
                                    ))
                         )
                       ),
                       column(
                         6,
                         box(
                           width = 12,
                           status = "danger",
                           solidHeader = TRUE,
                           DT::dataTableOutput("Modules") %>% withSpinner(color = "#FF0000")
                         )
                       )), 
                    
                    
                       )
                       ),
             
             
              tabPanel("Table Outputs",
                       fluidPage(
                         br(),
                       box(
                         width = 12,
                         status = "danger",
                         solidHeader = TRUE,
                       DT::dataTableOutput("Multi") %>% withSpinner(color = "#FF0000")
                       )
                       )),
              
              tabPanel("Summary Plots",
                       fluidPage(
                       br(),
                       
                      
                      # textOutput("clickevent"),
                       
                       box(title="Top Terms for each Module",width=12,status="danger",solidHeader = TRUE,
                       plotlyOutput("Bubble",width = "90%") %>% withSpinner(color = "#FF0000"),
                       uiOutput("bubble_databases"),
                       br(),
                       DTOutput("bubble_table")
                       ),
                      br(),
                           
                      # DTOutput("bubbleDT"),
                      
                       
                      box(title="Top Terms for each Database",width=12,status="danger",solidHeader = TRUE,
                       uiOutput("ontology_mod_ui"),
                       plotOutput("Ontology",width="80%") %>% withSpinner(color = "#FF0000")
                      )
                       )
                      ),
                      
                      tabPanel("Enrichment Heatmap",
                               fluidPage(
                                 br(),

                     dropdownButton(status = "danger",icon=icon("gears","font-awesome"),circle = FALSE,size="sm",label = "Modify heatmap",
                                    box(status = "danger",solidHeader = TRUE,
                           sliderInput("height","Height",min=0,max=2000,value=0),
                           sliderInput("width","Width",min=0,max=2000,value=0),
                           selectInput(inputId="color",
                                       label="Change heatmap color palette",
                                       choices=c("Default",COLORS),
                                       selected = "Default"
                           ),
                           sliderInput(inputId="breaks",
                                       label="Adjust Color Breaks",
                                       min=3,max=7,value=3,step = ),
                           prettyCheckbox(
                             inputId = "show_anno",
                             label = "Show Database Annotation", 
                             value = TRUE,
                             status = "danger",
                             shape = "curve"
                           )
                           
                           )),
                     br(),
                          plotOutput("heatmap") %>% withSpinner(color = "#FF0000")
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
      updateMultiInput(session,"select_mods",choices=unique(modules()$Modules))
    }
    
    if(input$run<1){
      hideTab("Tabs","Table Outputs")
    }else{showTab("Tabs","Table Outputs")}
    
    if(input$run<1){
      hideTab("Tabs","Summary Plots")
    }else{showTab("Tabs","Summary Plots")}
    
    if(input$run<1){
      hideTab("Tabs","Enrichment Heatmap")
    }else{showTab("Tabs","Enrichment Heatmap")}
    
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
    modules() %>% setNames(.,c("Gene Name","Module ID"))
    },filter="top",rownames=FALSE,options=list(pageLength=10)
  )



output$total_modules <- renderValueBox({
  all <- length(unique(modules()$Modules))
  active <- all-length(input$select_mods)
  
  valueBox(value =paste(active,"of",all),
           subtitle = "Active Modules",
           color = "red",
           icon = icon("network-wired",lib = "font-awesome")
           )
})

output$total_genes <- renderInfoBox({
  infoBox(value = length(unique(modules()$id)),
           title = "Unique Gene Names",
           color = "red",
           icon = icon("dna",lib = "font-awesome")
  )
})

# Module summary

output$modsummary_plot <-
 renderPlot({
  sum <- modules() %>%
     count(Modules)
  
  ggplot(sum, aes(n))+
     geom_histogram(fill="black", color="red",binwidth = 1)+
     xlab("Module Size")+
     ylab("Number of Modules")+
     theme_classic()

 })


  
  # Single DB Function
  
  Out <- eventReactive(input$do,{
   
        
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
         
        Output %>% mutate(SIG=-log10(P.value)) %>% arrange(SIG)
        }) 
  
  
		
   # Multiple DB Function
   
  Multi <- eventReactive(input$run,{

           modules <- read.csv(input$file$datapath,header=T)
           
           all.modules <- unique(modules$Modules)
           
 
           lapply(all.modules,function(x){
              genes <- modules %>% filter(Modules==x) %>% dplyr::pull(id)
              
           lapply(input$database,function(db){
                 output <- enrichr(genes,db)
                 Temp <- output[[1]] %>% arrange(Adjusted.P.value) # %>% .[1:input$rank,]
                 if(nrow(Temp)>0){
                 Temp %>% mutate(Database=db)
                 }
                 }) %>% bind_rows() %>% mutate(Module_ID=x) %>% mutate(Module.size=length(genes))

           }) %>% bind_rows() %>% select(Module_ID,
                                         Module.size,
                                         Database,
                                         Term,
                                         P.value,
                                         Adjusted.P.value,
                                         Odds.Ratio,
                                         Genes)
           

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
                          xaxis = list(showline= T, linewidth=2, linecolor='black'),
                          yaxis = list(showline= T, linewidth=2, linecolor='black'))
    
    fig
    
  })
   
  # clicked_module
  clicked_module <- reactiveValues(ID="")
   
  observeEvent(event_data("plotly_click"),{
    output$bubble_table <- renderDT({
   
    cords <- event_data("plotly_click")
    
    curve <- sort(input$database)[cords$curveNumber+1] 
    
      
    mod <- bubble_data() %>% {if(input$bubble_db=="all") dplyr::filter(.,Database==curve) else .} %>% .[cords$pointNumber+1,] %>% pull(Module_ID)
  
    clicked_module$ID <- mod
    #  mod <- Multi() %>%
    #    filter(round(Odds.Ratio,digits = 1)==cords$x) %>%
    #    pull(Module_ID)
      
   datatable(
     Multi() %>%
        filter(Module_ID==mod) %>%
        arrange(Adjusted.P.value) %>%
        select("Module ID"=Module_ID,
               Database,
               Term,
               "P-Value"=P.value,
               "Adjusted P-Value"=Adjusted.P.value,
               "Odds Ratio"=Odds.Ratio,Genes),
     rownames=FALSE,options=list(scrollX = TRUE),filter="top"
   ) %>% formatRound(c(4:6),digits=3) 
    })
    
    
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
  
  output$ontology_mod_ui <-renderUI({
    
    selectInput("ontology_mod","Select Module",choices=unique(Multi()$Module_ID),selected=unique(Multi()$Module_ID)[1])
  })
  
 
  
  output$Ontology <- renderPlot({
    req(is.null(input$ontol))
    
   data <- Multi() %>% filter(Module_ID==input$ontology_mod) %>% group_by(Database) %>%
     dplyr::filter(Adjusted.P.value==min(Adjusted.P.value)) %>% distinct(Module_ID,Adjusted.P.value,.keep_all=TRUE)
    
    ggplot(data,aes(x=-log10(Adjusted.P.value),y=reorder(Term,Adjusted.P.value)))+
      geom_bar(stat="identity",aes(fill=Database))+
     # geom_text(aes(label=Module_ID), position=position_dodge(width=0.9), vjust=-0.25)+
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
   
  
  
  
  output$heatmap <- renderPlot({
    
    Terms <- Multi() %>% group_by(Module_ID,Database) %>% filter(Adjusted.P.value==min(Adjusted.P.value) & Adjusted.P.value<0.05) %>% distinct(Module_ID,Adjusted.P.value,.keep_all=TRUE) %>%
      pull(Term)
    data <- Multi() %>% dplyr::filter(Term %in% Terms) %>%
      select(Term,Database,Module_ID,Adjusted.P.value) %>% mutate(Adjusted.P.value=-log10(Adjusted.P.value)) %>%
      pivot_wider(names_from =Module_ID,
                  values_from = Adjusted.P.value) %>% tibble::column_to_rownames(var = "Term")
    
    annotation <- data %>% select(Database)

    pheatmap(data %>% select(-Database) %>% replace(is.na(.), 0),
             border_color = "black",
             annotation_row = if(input$show_anno==TRUE){annotation}else{NA},
             color=col.pal()
             )

  },width=function(){if(input$width>0){input$width}else{"auto"}},height=function(){if(input$height>0){input$height}else{"auto"}})
  
  
  
  # For colors
  col.pal <- reactive({if(input$color=="Default"){colorRampPalette(rev(brewer.pal(n = input$breaks, name =
                                                                                    "RdYlBu")))(100)}else{colorRampPalette(rev(brewer.pal(n=input$breaks
                                                                                                                                          ,name= input$color)))(100)}
  })
  
  observeEvent(input$color,{if(input$color=="Default"){
    updateSliderInput(session, "breaks", label = "Number of color breaks in heatmap",
                      min=3,
                      max=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% "RdYlBu") 
                      %>% select(maxcolors) %>% .$maxcolors,
                      value=7
    )}
    else{
      
      
      
      
      updateSliderInput(session, "breaks", label = "Number of color breaks in heatmap",
                        min=3,
                        max=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% input$color) 
                        %>% select(maxcolors) %>% .$maxcolors,
                        value=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% input$color) 
                        %>% select(maxcolors) %>% .$maxcolors
      )}
  }) 
   
  # Reactive output for multi-db function
  observeEvent(input$run,{
    
    output$Multi <- DT::renderDT({datatable(Multi() %>% arrange(Adjusted.P.value) %>%
                                              select("Module ID"=Module_ID,
                                                     Database,
                                                     Term,
                                                     "P-Value"=P.value,
                                                     "Adjusted P-Value"=Adjusted.P.value,
                                                     "Odds Ratio"=Odds.Ratio,
                                                     Genes),
                                            rownames=FALSE,
                                            filter="top",
                                            options=list(pageLength=10)
                                            ) %>% formatRound(4:6,digits=3)
      }) 
    
    updateTabsetPanel(session, inputId="Tabs", selected="Summary Plots")
    
    
    shinyjs::show("DOWNLOAD_BUTTON")
  
  })
  
  observeEvent(input$file,ignoreInit = TRUE,{
  updateTabsetPanel(session, inputId="Tabs", selected="Your Modules")
  })
  
# Downloadable csvs of selected dataset ----
  
  # # Single DB output
  # output$downloadData <- 
  # 
  #   downloadHandler(
  #     filename = function() {
  #     paste("report", ".csv", sep = "")
  #   },
  #     content = function(file) {
  #     write.csv(Multi(), file, row.names = TRUE)
  #   }
  #   )
  
  # Multi DB output
	output$downloadData.two <- 
	  
	  downloadHandler(
      filename = function() {
      paste("RichieEnrichr_Results_",Sys.Date(), ".csv", sep = "")
    },
      content = function(file) {
      write.csv(Multi(), file, row.names = FALSE)
      }
    )

  
  
 }
  

shinyApp(ui, server)