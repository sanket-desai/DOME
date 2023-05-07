#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinyFiles)
library(dplyr)
library(data.table)

############## Added CSS effect to print errors in red color #############
appCSS <-
    ".mandatory_star { color: red; }
   #error { color: red; }"

############### To check mandatory Fields #################
fieldsMandatory <- c("Uniprot_Gene_Name", "Ref_AA", "Position")
 fieldsMandatory <- c("Uniprot_Gene_Name", "Position")
 
 
 labelMandatory <- function(label) {
     tagList(
         label,
         span("*", class = "mandatory_star")
     )
 }

fieldsMandatory1 <- c("filename")

labelMandatory1 <- function(label) {
    tagList(
        label,
        span("*", class = "mandatory_star")
    )
}

fieldsMandatory2 <- c("Uniprot_Gene_Name2")

labelMandatory2 <- function(label) {
    tagList(
        label,
        span("*", class = "mandatory_star")
    )
}

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cosmo"),
                shinyjs::useShinyjs(),
                shinyjs::inlineCSS(appCSS),
                navbarPage(
                    "DOME (DOmain Mutation Estimator)",
                    tabPanel("File Upload",
                             sidebarPanel(
                                 div(
                                     tags$h3(strong("Input:")),
                                     tags$br(),
                                     id = "form1",
                                     fileInput("filename", labelMandatory1("Choose your input file")),
                                     #textInput("Uniprot_Gene_Name", labelMandatory("Uniprot Gene Name"), "", placeholder = "e.g. EGFR_HUMAN"),
                                     #textInput("Ref_AA", labelMandatory("Reference Amino Acid"), "", placeholder = "e.g. L"),
                                     #textInput("Position", labelMandatory("Position of Amino Acid"), "", placeholder = "e.g. 858"),
                                     #textInput("Alt_AA", label = "Altered Amino Acid", "", placeholder = "e.g. R"),
                                     tags$br(),
                                     actionButton("submit1", "Submit", class = "btn-primary")
                                 ),
                                 shinyjs::hidden(
                                     div(
                                         id = "thankyou_msg1",
                                         h3("Thanks, your response was submitted successfully!"),
                                         actionLink("submit_another1", "Submit another Job")
                                     )
                                 ),
                                 shinyjs::hidden(
                                     span(id = "submit_msg1", "Submitting..."),
                                     div(id = "error1",
                                         div(br(), tags$b("Error: "), span(id = "error_msg1"))
                                     )
                                 )
                                 
                             ), ############## end of sidebarPanel
                             mainPanel(
                                 h1(strong("OUTPUT:")),
                                 br(),
                                 #h3("Output :"),
                                 dataTableOutput("domain_table1"),
                                 downloadButton("downloadData", "Save to disk")
                                 #downloadButton("download", "Download .tsv")
                                 #shinyjs::toggleState(id = "download", condition = mandatoryFilled)
                                 
                                 
                             )
                    )
                    
                    
                ))

# Define server logic required to draw a histogram
server <- function(input, output) {
       ############################### Navbar 2
    
    observe({
        mandatoryFilled1 <-
            vapply(fieldsMandatory1,
                   function(x) {
                       !any(is.null(input[[x]])) && input[[x]] != ""
                   },
                   logical(1))
        mandatoryFilled1 <- all(mandatoryFilled1)
        
        shinyjs::toggleState(id = "submit1", condition = mandatoryFilled1)
    })
    
    observeEvent(input$submit_another1, {
        shinyjs::reset("form1")
        shinyjs::show("form1")
        shinyjs::hide("thankyou_msg1")
    })
    
    observeEvent(input$submit1, {
        shinyjs::disable("submit1")
        shinyjs::show("submit_msg1")
        shinyjs::hide("error1")
        
        tryCatch({
            #saveData(formData())
            print(input$filename$name)
            file1 <- input$filename
            file1$filename$name
            #ugene_name <- input$Uniprot_Gene_Name
            #refaa <- input$Ref_AA
            #pos <- input$Position
            #altaa <- input$Alt_AA
            #print(ugene_name)
            #print(refaa)
            #print(pos)
            #print(altaa)
            domain1 <- NULL
            command_filtration1 <- paste("python scripts/metaaccessor.py" , input$filename$name, " > output/output_loop.csv", sep =" ")
            print(command_filtration1)
            system(command_filtration1, intern = TRUE)
            domain1 <- read.csv(file = 'output/output_loop.csv', header = T, sep = "\t")
            #print(domain1)
            domain1_1 <-subset(domain1,domain1$tier == 1)
            domain1_2 <-subset(domain1,domain1$tier == 2)
            domain1_3 <-subset(domain1,domain1$tier == 3)
            domain1_4 <-subset(domain1,domain1$tier == 4)
            domain1_0 <-subset(domain1,domain1$tier == 0)
            domain1_sort <- rbind(domain1_1,domain1_2,domain1_3,domain1_4,domain1_0)
            #domain1_sort <- domain1_sort[ -c(1,11,12,16,17,18,20,25,26,28,29) ]
            #domain1_sort <- domain1_sort[ -c(1,11,12,16,17,18,20,25,26,28,29) ]
            #print(domain1_sort)
            domain1_sort_filter <- mutate(domain1_sort, Status = ifelse((domain1_sort$tier == 1) & (domain1_sort$statscore == "H"), "Hotspot",ifelse((domain1_sort$tier == 1) & (domain1_sort$statscore == "R"), "Resistant",ifelse((domain1_sort$tier == 1) & (domain1_sort$statscore == "H,R"), "Hotspot / Resistant",ifelse((domain1_sort$meanscore >= 0.7), "High confidence",ifelse((domain1_sort$statscore > 0) & (domain1_sort$statscore < 0.7), "Moderate score", "-"))))))
            domain1_sort_filter_sorted <- domain1_sort_filter[,c(ncol(domain1_sort_filter),1:(ncol(domain1_sort_filter)-1))]
            domain1_sort_filter_sorted_clean <- domain1_sort_filter_sorted[ -c(2,12,13,17,18,19,21,26,27,29,30) ]
            #print(d)
            #d <- mutate(domain1_sort, uniprot = ifelse((domain1_sort$"uniprotsite" != "-") & (domain1_sort$"X3d_uniprotsite_close" == "-"), domain1_sort$"uniprotsite",ifelse((domain1_sort$"uniprotsite" == "-") & (domain1_sort$"X3d_uniprotsite_close" != "-"), domain1_sort$"X3d_uniprotsite_close",ifelse((domain1_sort$"phosphosite" != "-") & (domain1_sort$"X3d_pspsite_close" == "-"), domain1_sort$"phosphosite",ifelse((domain1_sort$"phosphosite" == "-") & (domain1_sort$"X3d_pspsite_close" != "-"), domain1_sort$"X3d_pspsite_close","-")))))
            #d <- mutate(d, hotspot = ifelse((d$"hotspot3dann" != "-") & (d$"X3d_hotspot_close" == "-"), d$"hotspot3dann",ifelse((domain1_sort$"tier" == 1) & (domain1_sort$"statscore" == "H"), "Hotspot")))
            #d <- mutate(d, analog = ifelse((d$"X3d_statsig_close" != "-") & (d$"X3d_analog_close" == "-"), d$"X3d_statsig_close",ifelse((d$"X3d_statsig_close" == "-") & (d$"X3d_analog_close" != "-"), d$"X3d_analog_close","-")))
            #domain_final <- d[-c(12,13,14,15,16,17,18,19)]
            colnames(domain1_sort_filter_sorted_clean) <- c("DOME Pred","Entry name","Position","Reference","Altered","Domain","Mutation type","Analogous to","CADD","COSMIC count","Clinvar","Site","Phosphosite","Proximity to significant","Proximity to functional site", "Proximity to PSP site", "Proximity to 3D hotspot","Is interface","Entropy","DOME Score")
            #output$domain_table1 <- renderDataTable(expr = domain1, options = list(scrollX=T)
            output$domain_table1 <- renderDataTable(expr = domain1_sort_filter_sorted_clean, options = list(scrollX=T)
            )
            # Downloadable csv of selected dataset ----
            observe({
                shinyjs::toggleState(id = "downloadData", !is.null(domain1_sort_filter_sorted_clean))
            })
            output$downloadData <- downloadHandler(
                filename = function() {
                    paste("output", ".csv", sep = "")
                },
                content = function(file) {
                    write.csv(domain1_sort_filter_sorted_clean, file, row.names = FALSE)
                }
            )
            Sys.sleep(2)
            
            
            #df <- as.data.frame(do.call(cbind, list(filname,variant_count)))
            #colnames(df) <- c("Input Files", "Varaints Counts")
            #print(df)
            #output$variants_count <- renderDataTable(expr = df)
            shinyjs::reset("form1")
            shinyjs::hide("form1")
            shinyjs::show("thankyou_msg1")
        },
        error = function(err) {
            shinyjs::html("error_msg1", err$message)
            shinyjs::show(id = "error1", anim = TRUE, animType = "fade")
        },
        finally = {
            shinyjs::enable("submit1")
            shinyjs::hide("submit_msg1")
        })
    })    
}

# Run the application
shinyApp(ui = ui, server = server)
