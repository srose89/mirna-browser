## Shiny miRNA Browser
## 3/25/2017

library(shiny)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # set title of page

  navbarPage(strong("ImmGen miRNA Browser", color = "red"), 
             tabPanel("miRNA Expression",
                      fluidRow(
                        column(12,
                               h3("Cell Type Expression", align = "center"),
                               plotOutput("exp_plot")
                        ),
                        fluidRow(theme = "paper.css",
                        column(3, class = "well",
                          h3("Cell Type", align = "left"), 
                          checkboxGroupInput("cell_type", label = "",
                                             choices = unique(names(sample_order)), 
                                             selected = "")),
                        
                        column(3, class = "well",
                               h3("Cell Class", align = "left"), 
                               checkboxGroupInput("cell_class", label = "",
                                                  choices = distinct(exp.lm, cell_class) %>%
                                                    collect() %>% .$cell_class, 
                                                  selected = "")),

                            column(3, class = "well",
                              h3("miRNA", align = "left"),
                              selectInput("miRNA", label = "", 
                                      choices = distinct(exp.lm, miRNA) %>% 
                                        arrange(as.character(miRNA)) %>% collect() %>% .$miRNA,
                                      selected = 'miR-223-3p')
                            ),
                        
                            column(3, class = "well", 
                              h3("Scale", align = "left"),
                              selectInput("exp", label = "", choices = c('Linear', 'Log2'), 
                                      selected = 'Log2')
                            )
                          )
                      )
              )

             ,tabPanel("Transcript-View",
                      # plot output section
                        fluidRow(
                          column(12,
                            h3("miRNA Regulation of Transcript", align = "center"),
                            br(),
                            plotOutput("tscan_plot", height = "500px"),
                            br()
                          )
                        ),
                        fluidRow(
                          column(4, class = "well", 
                          h4("Gene"),
                          textInput("gene", label = "", value = "",
                                      placeholder = "Ex. Myb"
                                      ),
                          uiOutput("tscript_display")
                          ),
                          # cell type choices
                          column(4, class = "well", 
                          h4("Cell Type"),
                          selectInput("cell_type.ts", label = "", 
                                      choices = distinct(exp.lm, cell_type) %>%
                                        select(cell_type) %>% collect(), selected = "B1aB"),
                          
                          selectInput("cell_type.ts2", label = "", 
                                      choices = distinct(exp.lm, cell_type) %>%
                                        select(cell_type) %>% collect(), selected = "B1aB")
                        ),
                        # create a column to select filtering and 2 panel view
                        column(4, class = "well",
                               h4("Parameters"),
                               sliderInput("exp_threshold", 
                                           h6("Expression threshold"), 
                                           min = 0, max = 15, value = 0),
                               sliderInput("context_threshold", 
                                           h6("context++ score threshold"), 
                                           min = 0, max = .2, value = 0),
                               checkboxInput("hd", strong("highlight differences"), 
                                             value = F),
                               sliderInput("diff_threshold", 
                                           h6("Expression difference threshold"), 
                                           min = 1, max = 10, value = 3, step = .5)
                        )
                        ), 
                       # download and target output
                       fluidRow(
                         column(12, 
                                h3("Target List"),
                                htmlOutput("tscan_text"),
                                br(), 
                                uiOutput("download_display"),
                                br(),
                                DT::dataTableOutput("tscan_list")
                                )
                       )
                      )#,
             # tabPanel("Help", 
             #          # this will be a page just explaining what is going on
             #          fluidRow(
             #            column(12, 
             #                   h3("Help Information"), 
             #                   br()#,
             #                   )
             #          ))
                    
  )
   , 
  title = "ImmGen miRNA Browser",
  theme = "paper.css")
)