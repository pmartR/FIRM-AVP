#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source("./functions/predict_AVP_sequences.R")
library(shiny)
library(markdown)
library(dplyr)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("FIRM-AVP: A Tool for Antiviral Peptide Prediction"),
    hr(),
    sidebarLayout(
        sidebarPanel(
            
            fileInput("fasta_file", "Choose FASTA file for Sequence Prediction",
                      multiple = FALSE
            ),
            textInput("fasta_text", "Enter a Sequence for Prediction", value = "", width = NULL, placeholder = NULL),
            
            fileInput("avp_fasta_file", "Add additional AVP Sequences to Training (FASTA)",
                      multiple = FALSE
            ),
            fileInput("nonavp_fasta_file", "Add additional Non-AVP Sequences to Training (FASTA)",
                      multiple = FALSE
            ),
            actionButton(inputId="go", label="Predict", icon = icon("calculator")),
            downloadButton("downloadResults", "Download Results")
        ),
        mainPanel(
            tabsetPanel(id = "inTabset",
                tabPanel("Welcome", id = "t1",
                         includeMarkdown("README.md"),
                         includeMarkdown("FASTA Formatting.md"),
                         downloadButton("downloadFasta", "Download Example")),
                #tabPanel("Predicted AVP Sequences",value = "t2", br(), htmlOutput("results"))
                tabPanel("Predicted AVP Sequences",value = "t2", br(),h3("Output Probabilities"), DT::dataTableOutput("results"))
                
            ) 
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$go, {
        updateTabsetPanel(session,"inTabset",selected ="t2")
    })
    
    # allow an example FASTA for download
    output$downloadFasta <- downloadHandler(
        filename = function() {
            paste("Example", ".fasta", sep = "")
        },
        content = function(file) {
            lines = readLines('./data/input_seq.fasta')
            write.table(lines, file, row.names = FALSE)
        }
    )
    # if there are new sequences to predict, add them, otherwise use default
    testing_data <- reactive({
        if(is.null(input$fasta_file)){
            if(input$fasta_text == ""){
                testing_file = read.csv("./data/input_seq.csv", header=TRUE)
                return(testing_file)
            } else {
                # try and parse the input text
                return(feature_extraction(fasta_file_path=NA, aa_idx_path = "./data/AAidx.csv", raw_sequences=toupper(trimws(input$fasta_text))))
            }
        }
        else{
            return(feature_extraction(input$fasta_file$datapath))
        }
    })
    
    training_data <- reactive({
        #old_sequences <- read.csv("./data/selected_train_test_merged_file.csv", header = TRUE)
        old_sequences <- read.csv("./data/non_sars_merged_train_test.csv", header = TRUE)
        
        features_selected=c("aac_1","aac_2","aac_3","aac_4","aac_6","aac_7","aac_8","aac_9","aac_10","aac_11","aac_12","aac_15","aac_16","aac_17","aac_18","aac_19","aac_20","dipep_32","dipep_51","dipep_111","dipep_211","dipep_220","dipep_340","pseudo_1","pseudo_2","pseudo_3","pseudo_4","pseudo_5","pseudo_10","pseudo_11","pseudo_12","pseudo_14","pseudo_16","pseudo_18","pseudo_20","pseudo_21","pseudo_22","pseudo_23","pseudo_24","pseudo_25","amphipseudo_21","amphipseudo_22","amphipseudo_23","amphipseudo_24","amphipseudo_25","amphipseudo_26","amphipseudo_27","amphipseudo_29","amphipseudo_30","comp_1","comp_2","comp_3","comp_4","comp_5","comp_6","comp_10","comp_11","comp_13","comp_14","comp_15","comp_16","comp_17","comp_18","comp_19","comp_21","comp_22","comp_23","comp_24","tran_1","tran_2","tran_3","tran_4","tran_5","tran_6","tran_11","tran_12","tran_13","tran_14","tran_16","tran_17","tran_18","tran_19","tran_20","tran_21","tran_22","tran_23","tran_24","dist_1","dist_2","dist_3","dist_4","dist_7","dist_8","dist_9","dist_10","dist_11","dist_12","dist_13","dist_14","dist_15","dist_16","dist_17","dist_18","dist_22","dist_23","dist_24","dist_25","dist_26","dist_27","dist_28","dist_29","dist_30","dist_32","dist_34","dist_38","dist_41","dist_46","dist_47","dist_50","dist_52","dist_53","dist_55","dist_56","dist_61","dist_62","dist_63","dist_65","dist_67","dist_68","dist_70","dist_71","dist_72","dist_73","dist_76","dist_77","dist_78","dist_79","dist_82","dist_83","dist_84","dist_85","dist_86","dist_87","dist_88","dist_89","dist_90","dist_91","dist_93","dist_94","dist_97","dist_99","dist_100","dist_102","dist_103","dist_105","dist_106","dist_107","dist_108","dist_109","dist_112","dist_113","dist_114","dist_115","dist_116","dist_117","dist_118","dist_119","dist_120","ss_1","Output")
        old_sequences <- old_sequences[, which(colnames(old_sequences) %in% features_selected)]
        new_sequences <- NULL
        if(!is.null(input$avp_fast_file)){
            new_sequences  <- add_new_sequences(input$fasta_file$datapath, outcome_value = 1)
        }
        if(!is.null(input$nonavp_fasta_file)){
            new_sequences_nonavp  <- add_new_sequences(input$fasta_file$datapath, outcome_value = -1)
            if(!is.null(new_sequences)){
                rbind(new_sequences, new_sequences_nonavp)
            }else{
                new_sequences <- new_sequences_nonavp
            }
        }
        if(is.null(new_sequences)){
            return(old_sequences)
        }else{
            training_file <- rbind(new_sequences, old_sequences)
            return(training_file)
        }
        
    })
    
    results <- eventReactive(input$go,{
        if(is.null(input$fasta_file$datapath)){
            if(input$fasta_text == ""){
                withProgress( expr = {temp <- predict_AVP_sequences(training_data(),
                                                                    testing_data(),
                                                                    fasta_file_path = "./data/input_seq.fasta")
                incProgress(1)
                }, message = "Calculating...please wait")  
            } else {
                withProgress( expr = {temp <- predict_AVP_sequences(training_data(),
                                                                    testing_data(),
                                                                    fasta_file_path = NA,
                                                                    raw_sequence = toupper(trimws(input$fasta_text)))
                incProgress(1)
                }, message = "Calculating...please wait") 
            }
            
            #return(paste("Example file predictions: ", temp, sep = ""))
            return(temp)
        } else if(!is.null(input$fasta_file$datapath)){
            
            withProgress( expr = {temp <- predict_AVP_sequences(training_data(),
                                                                testing_data(),
                                                                fasta_file_path = input$fasta_file$datapath)
            incProgress(1)
            }, message = "Calculating...please wait")
            return(temp)
        }
    })
    
    #output$results <- renderUI({HTML(results())})

    output$results <- DT::renderDataTable(DT::datatable(results(), rownames = FALSE) %>% DT::formatStyle(
        'AVP',
        target = 'row',
        backgroundColor = DT::styleInterval(c(0.5), c('white', '#99CCFF'))
    )%>%
        DT::formatRound('AVP', 4) %>%
        DT::formatRound('Non-AVP', 4)
    )
    
    output$downloadResults <- downloadHandler(
        filename = function() {
            paste("AVPR_Results", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(results(), file, row.names = FALSE)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
