## server.R
## 3/25/2017


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  # output for standard miRNA expression barplot
  output$exp_plot <- renderPlot({
    exp <- switch(input$exp,
                  "Linear" = exp.lm,
                  "Log2" = exp.l2m)

    mirnaBarplot(exp, input$miRNA, input$cell_type, input$cell_class)
  })


  # enable choice of transcript ID if there is more than one for
  # a given gene

  tscript_choice <- reactive({
    req(input$gene)
    g <- firstup(input$gene)
    # added a req in case the gene name doesn't match
    num_tscript <- filter(utr, Gene_symbol == g & `Representative_transcript?` == 1) %>%
      collect() %>% .$Transcript_ID %>% length() %>% req()
    if(num_tscript > 1)
      return(TRUE)

    return(FALSE)
  })

  output$tscript_display <- renderUI({
    req(tscript_choice())
    g <- firstup(input$gene)
    # added a req in case the gene name doesn't match
    rep_tscript <- filter(utr, Gene_symbol == g & `Representative_transcript?` == 1) %>% collect() %>% .$Transcript_ID %>% req()
    # if there isn't just one representative transcript select the appropriate one
    selectInput("tscript_selection", "Select Transcript",
                  choices = rep_tscript, multiple = F)
  })

  # function to get transcript that will be used that I can then use for other things
  tscript_use <- reactive({
    #req(tscript_choice() == F | input$tscript_selection)
    g <- firstup(input$gene)
    #print("getting to tscript use")
    if(tscript_choice()){
      return(input$tscript_selection)
    }
    else{
      return(filter(utr, Gene_symbol == g & `Representative_transcript?` == 1) %>%
               collect() %>% .$Transcript_ID)
    }
  })

  # output for transcript level plot
  output$tscan_plot <- renderPlot({

    req(tscanTable())
    tscan_bar()

  })

  # generate plot for UTR for download or on screen output
  tscan_bar <- reactive({
    req(tscanTable())


    p <- plotGeneTscan(tscript_use(), tscanTable(),  utr,
                       exp_thresh = input$exp_threshold,
                       exp_diff_check = input$hd & (input$cell_type.ts != input$cell_type.ts2),
                       exp_diff_thresh = input$diff_threshold,
                       cell_types_chosen = c(input$cell_type.ts2,input$cell_type.ts))

    p
  })


  ## generate a table based on the gene and cell type input that will be used
  ## for display and downloading
  ## adding in the second cell type
  tscanTable <- reactive({

    # require that the representative transcript has been picked
    req(tscript_use())
    rep_tscript <- tscript_use()

    # filter expression table
    exp.f <- filter(exp.l2m, cell_type %in% c(local(input$cell_type.ts),
                                              local(input$cell_type.ts2))) %>%
      group_by(tscan_id, cell_type) %>%
      summarise(mean_expression = mean(expression))  %>% collect()
    # filter targetscan file
    ts <- collect(filter(tscan.full, Transcript_ID == rep_tscript &
                           (abs(`context++_score`) > local(input$context_threshold) |
                              is.na(`context++_score`))) %>%
      dplyr::select(miRNA, Site_Type, UTR_start, UTR_end, `context++_score`, `context++_score_percentile`))
    # join the two
    te.j <- left_join(ts, exp.f, by = c("miRNA" = "tscan_id"))
    # calculate differences if the box is checked
    if(input$hd == TRUE & (input$cell_type.ts != input$cell_type.ts2)){
      # change the character strings to symbols so they can be used in the
      # dplyr function
      cell1 <- sym(input$cell_type.ts)
      cell2 <- sym(input$cell_type.ts2)
      # print(head(na.omit(te.j) %>%
      #              tidyr::spread(cell_type, mean_expression) %>%
      #              mutate(exp_diff = round((!!cell1) - (!!cell2), 3))))
      te.j <- na.omit(te.j) %>%
        tidyr::spread(cell_type, mean_expression) %>%
        mutate(exp_diff = round((!!cell1) - (!!cell2), 3)) %>%
        reshape2::melt(c(1:6,9), c(7,8), variable.name = "cell_type",
             value.name = "mean_expression") %>%
        left_join(te.j, ., by = c("miRNA", "Site_Type", "UTR_start", "UTR_end",
                                  "context++_score","context++_score_percentile",
                                  "cell_type", "mean_expression"))

    }
    # return the processed table

    te.j
  })

  # render the dataTable as an output
  output$tscan_list <- DT::renderDataTable({
    tscanTable()
  }, rownames = F)


  # framework to display nothing if no gene has been put in
  output$download_display <- renderUI({
    if(nrow(tscanTable()) == 0)
      return("")
    fluidRow(
    downloadButton("tscanDownload", label = "Download Target Table"),
    downloadButton("tscanPlot", label = "Download UTR Plot"),
    numericInput("plotWidth", "Plot width",  10, min = 1, max = 30, width = '10%' ),
    numericInput("plotHeight", "Plot height",  8, min = 1, max = 30, width = '10%')
    )
  })

  # output for targetscan gene info
  output$tscan_text <- renderText({
    g <- firstup(input$gene)
    tab <- collect(filter(utr, Gene_Symbol == g))
    HTML(paste("<b>Gene Symbol:</b>", tab$Gene_symbol[1],
               "\n<b>Transcript ID:</b>", tab$Transcript_ID[1]))
  })

  ## output for downloadable csv of target sites and expression
  output$tscanDownload <- downloadHandler(
    filename = function() {
      cells.dl <- as.character(paste(unique(as.character(c(input$cell_type.ts, input$cell_type.ts2))), collapse = "_"))
      paste(paste("target_sites", input$gene, cells.dl,
                  sep = "_"),
            ".csv",sep = "")
    },
    content = function(file){
      write.csv(tscanTable(), file,  row.names = F, quote = F)
    }
  )
  ## output for downloadable pdf of target sites over a UTR
  output$tscanPlot <- downloadHandler(
    filename = function() {
      cells.dl <- as.character(paste(unique(as.character(c(input$cell_type.ts, input$cell_type.ts2))), collapse = "_"))
      paste(paste("target_sites", input$gene, cells.dl,
                  sep = "_"),
            ".pdf",sep = "")

    },
    content = function(file){
      # function to output plot display to pdf file
      ggsave(file, tscan_bar(),
             width = input$plotWidth,
             height = input$plotHeight,
             device = "pdf")
    }
  )

})
