# Load necessary libraries
library(shiny)
library(DT) # for data table display
library(dplyr)
library(openxlsx)
library(readxl)
library(tools)
library(bslib)
library(jsonlite)

list_top_genes <- fromJSON("/Users/bli/Documents/Projects/WebApps/pdac_molecular_subtyping_shiny_app/list_top_genes.json")

source("/Users/bli/Documents/Projects/WebApps/pdac_molecular_subtyping_shiny_app/utility_functions.R")

# days_on_study_colname <- "Days On Study"
options(shiny.maxRequestSize = 100 * 1024 ^ 2)

############
library(shiny)
library(readxl)
library(DT)
library(bslib)
library(jsonlite)

# Define UI for application with Minty theme
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "minty"),  # Green color scheme
  
  titlePanel("PDAC Molecular Subtyping"),
  
  # Custom CSS for multi-column checkbox layout
  tags$head(
    tags$style(HTML("
            .multicol {
                -webkit-column-count: 3; /* Chrome, Safari, Opera */
                -moz-column-count: 3;    /* Firefox */
                column-count: 3;
                -webkit-column-gap: 20px;
                -moz-column-gap: 20px;
                column-gap: 20px;
            }
        "))
  ),
  
  # Tabset panel with two tabs: Data and Analysis
  tabsetPanel(
    tabPanel("Data",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file1", "Choose Normalized bulk RNA-Seq expression data file (*.xlsx)",  # Keep Excel upload functionality
                           accept = c(".xlsx")),
                 checkboxInput("select_all", "Unselect All Samples", value = FALSE),
                 div(
                   class = "multicol",  
                   uiOutput("sample_selector")
                 )
               ),
               mainPanel(
                 DTOutput("table")
               )
             )
    ),
    
    tabPanel("Classification",
             sidebarLayout(
               sidebarPanel(
                 # Multi-choice checkbox for classifiers
                 checkboxGroupInput("classifiers", "Choose classifiers:",
                                    # choices = c("Moffitt", "PurIST", "Collisson", "Bailey", "Puleo", "Chan-Seng-Yue")),
                                    choices = c("Moffitt", "PurIST", "Collisson", "Bailey", "Puleo")),
                 
                 # Dropdown for selecting one classifier from the selected ones
                 uiOutput("order_by_selector"),  # Dynamic UI for "Order classification by"
                 
                 actionButton("run_analysis", "Run Analysis")  # Button for running analysis, 
               ),
               mainPanel(
                 verbatimTextOutput("analysis_output"), 
                 downloadButton("download_results", "Download Predictive Results"), # add a download button 
                 downloadButton("download_genedata", "Download Gene Data") # add another download button for gene profile data
               )
             )
    ), 
    
    tabPanel("Plot", 
             sidebarLayout(
               sidebarPanel(
                 
                 # Optional input for showing only N most upregulated genes
                 numericInput("n_genes", 
                              "Show only top N signature genes of each subtype:", 
                              value = NA, 
                              min = 1, 
                              step = 1),
                 
                 actionButton("show_heatmap", "Generate Analysis Plot")  # Button for running analysis, 
               ),
               mainPanel(
                 plotOutput("heatmap_plot", height = "1200px", width = "900px")
               )
             )
             
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  
  # create reactive value obj to store variables
  rv <- reactiveValues()
  
  # Reactive function to read the uploaded Excel file
  data <- reactive({
    req(input$file1)
    
    # Read the Excel file
    df <- read_excel(input$file1$datapath)
    
    # Check if the first column is named 'Hugo_Symbol'
    if (colnames(df)[1] != "Hugo_Symbol") {
      showNotification("Error: The first column must be named 'Hugo_Symbol'.", type = "error")
      return(NULL)
    }
    
    return(df)
  })
  
  # Create dynamic UI for column (sample) selection, excluding 'Hugo_Symbol'
  output$sample_selector <- renderUI({
    req(data())
    cols <- colnames(data())[-1]
    checkboxGroupInput("selected_samples", 
                       "Select Samples (Columns):", 
                       choices = cols,
                       selected = if (input$select_all) character(0) else cols)
  })
  
  # Observe when the "select_all" checkbox changes to update the sample selector
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "selected_samples", 
                             selected = if (input$select_all) character(0) else colnames(data())[-1])
  })
  
  # Reactive function to rename 'Hugo_Symbol' to 'gene_name' and extract selected data for analysis
  selected_data <- reactive({
    req(data(), input$selected_samples)
    df <- data()
    colnames(df)[1] <- "gene_name"
    df[, c("gene_name", input$selected_samples), drop = FALSE]
  })
  
  # Render DataTable with 'gene_name' plus selected samples
  output$table <- renderDT({
    req(selected_data())
    datatable(selected_data(), options = list(pageLength = 10))
  })
  
  # Dynamic UI for selecting "Order classification by"
  output$order_by_selector <- renderUI({
    req(input$classifiers)  # Ensure classifiers have been selected
    selectInput("order_by", "Order classification by:", choices = input$classifiers)
  })
  
  # Perform analysis when the button is clicked
  observeEvent(input$run_analysis, {
    # Ensure that a classifier has been selected for ordering
    req(input$classifiers, input$order_by)
    
    # convert selected_data() to data.frame
    dat_exp <- as.data.frame(selected_data())
    
    ## construct dat_pred
    for (idx in 1:length(input$classifiers)) {
      
      name_cls <- input$classifiers[idx]
      
      if (name_cls %in% c("Moffitt", "PurIST")) {
        
        dat_tmp <- predict_pdac_subtype_ss(dat_exp, name_cls)
        
      } else {
        
        dat_tmp <- predict_pdac_subtype(dat_exp, list_top_genes[[name_cls]])
        
      }
      
      if (idx == 1) {
        
        dat_pred <- dat_tmp %>%
          select(sample_id, Prediction)
        colnames(dat_pred) <- c("sample_id", name_cls)
        
      } else {
        
        dat_tmp <- dat_tmp %>%
          select(Prediction)
        colnames(dat_tmp) <- c(name_cls)
        
        dat_pred <- cbind(dat_pred, dat_tmp)
        
      }
    }
    # rownames(dat_pred) <- NULL
    
    # order by input$order_by
    dat_pred_ordered <- dat_pred[order(dat_pred[[input$order_by]]), ]
    rownames(dat_pred_ordered) <- NULL
  
    rv$dat_pred <- dat_pred
    rv$dat_exp <- dat_exp
    rv$dat_pred_ordered <- dat_pred_ordered
     
    output$analysis_output <- renderPrint({
      # analysis_result
      
      dat_pred_ordered
      
    })
    
  })
  
  
  # Define the download handler for exporting prediction/classification results
  output$download_results <- downloadHandler(
    filename = function() {
      paste("classification_results_", Sys.Date(), "_", format(Sys.time(), "%H%M%S"), ".xlsx", sep = "")
    },
    content = function(file) {
      req(rv$dat_pred_ordered)  # Ensure the data exists
      
      # Create a workbook
      wb <- createWorkbook()
      
      # Add dat_pred_ordered as the first sheet
      addWorksheet(wb, "Summary")
      writeData(wb, "Summary", rv$dat_pred_ordered)
      
      # Iterate through classifiers and add dat_tmp as separate sheets
      for (idx in 1:length(input$classifiers)) {
        name_cls <- input$classifiers[idx]
        
        if (name_cls %in% c("Moffitt", "PurIST")) {
          dat_tmp <- predict_pdac_subtype_ss(rv$dat_exp, name_cls)
        } else {
          dat_tmp <- predict_pdac_subtype(rv$dat_exp, list_top_genes[[name_cls]])
        }
        
        addWorksheet(wb, name_cls)
        
        # replacing NA with "NA"
        dat_tmp[is.na(dat_tmp)] <- "NA"
        
        writeData(wb, name_cls, dat_tmp)
      }
      
      # Save the workbook to the specified file
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # Define download handler for exporting genedata results
  output$download_genedata <- downloadHandler(
    
    filename = function() {
      paste0("genedata_results_", Sys.Date(), "_", format(Sys.time(), "%H%M%S"), ".xlsx", sep = "")
    }, 
    content = function(file) {
      
      # create a workbook
      wb <- createWorkbook()
      
      # iterate through classifiers and to sheets
      for (idx in 1:length(input$classifiers)) {
        name_cls <- input$classifiers[idx]
        
        if (name_cls %in% c("Moffitt", "PurIST")) {
          genedat_tmp <- generate_pdac_subtype_genedata_ss(rv$dat_exp, name_cls)
        } else {
          genedat_tmp <- generate_pdac_subtype_genedata(rv$dat_exp, list_top_genes[[name_cls]])
        }
        
        addWorksheet(wb, name_cls)
        
        # replacing NA with "NA"
        genedat_tmp[is.na(genedat_tmp)] <- "NA"
        
        writeData(wb, name_cls, genedat_tmp)
      }
      
      # save workbook to specified file
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # generate heatmap plot when the button is clicked
  observeEvent(input$show_heatmap, {
    
    # 
    # print(length(rv$dat_exp$gene_name))
    # print(length(unique(rv$dat_exp$gene_name)))
    # 
    # 
    
    if (input$order_by %in% c("Moffitt", "PurIST")) { # use all genes of moffitt and purist regardless of input$n_genes assigned
      assigner_genes <- as.vector(unlist(list_top_genes[[input$order_by]]))
      # print(assigner_genes)
    } else {
      
      if (!is.na(input$n_genes)) {
        
        # assigner_genes <- list_top_genes[[input$order_by]][[1]]$gene_name[1:100]
        
        # select top N most upregulated genes of each subtype and remove duplicates
        tmp_genes <- do.call(c, lapply(names(list_top_genes[[input$order_by]]), function(cls) {
          tmp <- list_top_genes[[input$order_by]][[cls]]
          tmp2 <- tmp[order(-tmp$value), ][1:input$n_genes, ]
          tmp2$gene_name
        }))
        
        assigner_genes <- tmp_genes[!duplicated(tmp_genes)]
      
      } else {
        
       assigner_genes <- list_top_genes[[input$order_by]][[1]]$gene_name
      
      }
    }
    
    dat_exp_plot <- rv$dat_exp %>%
      filter(gene_name %in% assigner_genes)
    dat_exp_plot$gene_name <- factor(dat_exp_plot$gene_name, levels = assigner_genes)
      
    output$heatmap_plot <- renderPlot(
      generate_heatmap_with_annotations(
        dat_exp_plot,
        # rv$dat_exp, 
        rv$dat_pred_ordered, 
        gene_name_column_string = "gene_name", 
        order_by = input$order_by, 
        gene_limit = 100, 
        order_genes_by = assigner_genes
      )
    )
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)