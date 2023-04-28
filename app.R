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
#fieldsMandatory <- c("Uniprot_Gene_Name", "Ref_AA", "Position")
# fieldsMandatory <- c("Uniprot_Gene_Name", "Position")
# 
# 
# labelMandatory <- function(label) {
#     tagList(
#         label,
#         span("*", class = "mandatory_star")
#     )
# }

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
                    # tabPanel("Search Panel",
                    #          sidebarPanel(
                    #              div(
                    #                  tags$h3(strong("Input:")),
                    #                  tags$br(),
                    #                  id = "form",
                    #                  textInput("Uniprot_Gene_Name", labelMandatory("Uniprot Gene Name"), "", placeholder = "e.g. EGFR_HUMAN"),
                    #                  textInput("Ref_AA", label="Reference Amino Acid", "", placeholder = "e.g. L"),
                    #                  textInput("Position", labelMandatory("Position of Amino Acid"), "", placeholder = "e.g. 858"),
                    #                  textInput("Alt_AA", label = "Altered Amino Acid", "", placeholder = "e.g. R"),
                    #                  tags$br(),
                    #                  actionButton("submit", "Submit", class = "btn-primary")
                    #              ),
                    #              shinyjs::hidden(
                    #                  div(
                    #                      id = "thankyou_msg",
                    #                      h3("Thanks, your response was submitted successfully!"),
                    #                      actionLink("submit_another", "Submit another Job")
                    #                  )
                    #              ),
                    #              shinyjs::hidden(
                    #                  span(id = "submit_msg", "Submitting..."),
                    #                  div(id = "error",
                    #                      div(br(), tags$b("Error: "), span(id = "error_msg"))
                    #                  )
                    #              )
                    #              
                    #          ), ############## end of sidebarPanel
                    #          mainPanel(
                    #              h1(strong("OUTPUT:")),
                    #              br(),
                    #              #h3("Output :"),
                    #              dataTableOutput("domain_table")
                    #              #downloadButton("download", "Download .tsv")
                    #              #shinyjs::toggleState(id = "download", condition = mandatoryFilled)
                    #              
                    #              
                    #          )
                    # ),
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
                    # tabPanel("Search Gene",
                    #          sidebarPanel(
                    #              div(
                    #                  tags$h3(strong("Input:")),
                    #                  tags$br(),
                    #                  id = "form2",
                    #                  textInput("Uniprot_Gene_Name2", labelMandatory2("Uniprot Gene Name"), "", placeholder = "e.g. EGFR_HUMAN"),
                    #                  #textInput("Ref_AA2", label="Reference Amino Acid", "", placeholder = "e.g. L"),
                    #                  textInput("startPosition", label="Start Position of Amino Acid", "", placeholder = "e.g. 858"),
                    #                  textInput("endPosition", label="End Position of Amino Acid", "", placeholder = "e.g. 900"),
                    #                  #textInput("Alt_AA2", label = "Altered Amino Acid", "", placeholder = "e.g. R"),
                    #                  tags$br(),
                    #                  actionButton("submit2", "Submit", class = "btn-primary")
                    #              ),
                    #              shinyjs::hidden(
                    #                  div(
                    #                      id = "thankyou_msg2",
                    #                      h3("Thanks, your response was submitted successfully!"),
                    #                      actionLink("submit_another2", "Submit another Job")
                    #                  )
                    #              ),
                    #              shinyjs::hidden(
                    #                  span(id = "submit_msg2", "Submitting..."),
                    #                  div(id = "error2",
                    #                      div(br(), tags$b("Error: "), span(id = "error_msg2"))
                    #                  )
                    #              )
                    #              
                    #          ), ############## end of sidebarPanel
                    #          mainPanel(
                    #              h1(strong("OUTPUT:")),
                    #              br(),
                    #              #h3("Output :"),
                    #              dataTableOutput("domain_table2"),
                    #              downloadButton("downloadData2", "Save to disk")
                    #              #downloadButton("download", "Download .tsv")
                    #              #shinyjs::toggleState(id = "download", condition = mandatoryFilled)
                    #              
                    #              
                    #          )
                    # )
                    
                    
                ))

# Define server logic required to draw a histogram
server <- function(input, output) {
# /*
#         observe({
#         mandatoryFilled <-
#             vapply(fieldsMandatory,
#                    function(x) {
#                        !is.null(input[[x]]) && input[[x]] != ""
#                    },
#                    logical(1))
#         mandatoryFilled <- all(mandatoryFilled)
#         
#         shinyjs::toggleState(id = "submit", condition = mandatoryFilled)
#         shinyjs::toggleState(id = "downloadData", condition = mandatoryFilled)
#     })
#     
#     observeEvent(input$submit_another, {
#         shinyjs::reset("form")
#         shinyjs::show("form")
#         shinyjs::hide("thankyou_msg")
#     })
#     observeEvent(input$submit, {
#         shinyjs::disable("submit")
#         shinyjs::show("submit_msg")
#         shinyjs::hide("error")
#         
#         tryCatch({
#             #saveData(formData())
#             ugene_name <- input$Uniprot_Gene_Name
#             refaa <- input$Ref_AA
#             pos <- input$Position
#             altaa <- input$Alt_AA
#             print(ugene_name)
#             print(refaa)
#             print(pos)
#             print(altaa)
#             command_filtration <- paste("python scripts/metaaccessor.py" , paste(ugene_name, pos, altaa, sep = ","), sep =" ")
#             print(command_filtration)
#             system(command_filtration, intern = TRUE)
#             domain <- read.csv(file = 'output/output1.csv', header = T)
#             output$domain_table <- renderDataTable(expr = domain, options = list(scrollX=T)
#             #output$domain_table <- renderDataTable(expr = domain
#             )
#             Sys.sleep(2)
#             
#             
#             #df <- as.data.frame(do.call(cbind, list(filname,variant_count)))
#             #colnames(df) <- c("Input Files", "Varaints Counts")
#             #print(df)
#             #output$variants_count <- renderDataTable(expr = df)
#             shinyjs::reset("form")
#             shinyjs::hide("form")
#             shinyjs::show("thankyou_msg")
#         },
#         error = function(err) {
#             shinyjs::html("error_msg", err$message)
#             shinyjs::show(id = "error", anim = TRUE, animType = "fade")
#         },
#         finally = {
#             shinyjs::enable("submit")
#             shinyjs::hide("submit_msg")
#         })
#     })
 
       ############################### Navbar 2
    
    observe({
        mandatoryFilled1 <-
            vapply(fieldsMandatory1,
                   function(x) {
                       !is.null(input[[x]]) && input[[x]] != ""
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
            command_filtration1 <- paste("python scripts/metaaccessor.py" , input$filename$name, " >> output/output_loop.csv", sep =" ")
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
            domain1_sort <- domain1_sort[ -c(1,11,12,16,17,18,20,25,26,28,29) ]
            #d <- mutate(domain1_sort, uniprot = ifelse((domain1_sort$"uniprotsite" != "-") & (domain1_sort$"X3d_uniprotsite_close" == "-"), domain1_sort$"uniprotsite",ifelse((domain1_sort$"uniprotsite" == "-") & (domain1_sort$"X3d_uniprotsite_close" != "-"), domain1_sort$"X3d_uniprotsite_close",ifelse((domain1_sort$"phosphosite" != "-") & (domain1_sort$"X3d_pspsite_close" == "-"), domain1_sort$"phosphosite",ifelse((domain1_sort$"phosphosite" == "-") & (domain1_sort$"X3d_pspsite_close" != "-"), domain1_sort$"X3d_pspsite_close","-")))))
            #d <- mutate(d, hotspot = ifelse((d$"hotspot3dann" != "-") & (d$"X3d_hotspot_close" == "-"), d$"hotspot3dann",ifelse((d$"hotspot3dann" == "-") & (d$"X3d_hotspot_close" != "-"), d$"X3d_hotspot_close","-")))
            #d <- mutate(d, analog = ifelse((d$"X3d_statsig_close" != "-") & (d$"X3d_analog_close" == "-"), d$"X3d_statsig_close",ifelse((d$"X3d_statsig_close" == "-") & (d$"X3d_analog_close" != "-"), d$"X3d_analog_close","-")))
            #domain_final <- d[-c(12,13,14,15,16,17,18,19)]
            colnames(domain1_sort) <- c("Entry name","Position","Reference","Altered","Domain","Mutation type","Analogous to","CADD","COSMIC count","Clinvar","Site","Phosphosite","Proximity to significant","Proximity to functional site", "Proximity to PSP site", "Proximity to 3D hotspot","Is interface","Entropy","DOME Score")
            #output$domain_table1 <- renderDataTable(expr = domain1, options = list(scrollX=T)
            output$domain_table1 <- renderDataTable(expr = domain1_sort, options = list(scrollX=T)
            )
            # Downloadable csv of selected dataset ----
            observe({
                shinyjs::toggleState(id = "downloadData", !is.null(domain1_sort))
            })
            output$downloadData <- downloadHandler(
                filename = function() {
                    paste("output_rohit", ".csv", sep = "")
                },
                content = function(file) {
                    write.csv(domain1, file, row.names = FALSE)
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
    ############################### Navbar 3
    
    # observe({
    #     mandatoryFilled2 <-
    #         vapply(fieldsMandatory2,
    #                function(x) {
    #                    !is.null(input[[x]]) && input[[x]] != ""
    #                },
    #                logical(1))
    #     mandatoryFilled2 <- all(mandatoryFilled2)
    #     
    #     shinyjs::toggleState(id = "submit2", condition = mandatoryFilled2)
    #     shinyjs::toggleState(id = "downloadData2", condition = mandatoryFilled2)
    # })
    # 
    # observeEvent(input$submit_another2, {
    #     shinyjs::reset("form")
    #     shinyjs::show("form2")
    #     shinyjs::hide("thankyou_msg2")
    # })
    # 
    # observeEvent(input$submit2, {
    #     shinyjs::disable("submit2")
    #     shinyjs::show("submit_msg2")
    #     shinyjs::hide("error2")
    #     
    #     tryCatch({
    #         #saveData(formData())
    #         #print(input$filename$name)
    #         #file1 <- input$filename
    #         #file1$filename$name
    #         ugene_name2 <- input$Uniprot_Gene_Name2
    #         #refaa <- input$Ref_AA
    #         pos_start <- input$startPosition
    #         pos_end <- input$endPosition
    #         #altaa <- input$Alt_AA
    #         #print(ugene_name)
    #         #print(refaa)
    #         #print(pos)
    #         #print(altaa)
    #         domain2 <- NULL
    #         #print("Hello")
    #         #header <- data.frame("protein","mutpos","refaa","altaa","domain","statscore","analogtomut","caddscore","cosmiccount","clinvarsig","uniprotsite","phosphosite","3d_statsig_close","3d_uniprotsite_close","3d_pspsite_close","3d_hotspot_close","interfaceposition","entropyscore","mean4")
    #         #write.table(header,file="output/output.csv",quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
    #         #print(input$endPosition)
    #         if (input$endPosition==""){
    #           command_filtrationh <- paste("zcat data/domerefdata/final/proteomedomepred.tsv.gz | head -1 > output/output.csv")
    #           command_filtration2 <- paste("tabix data/domerefdata/final/proteomedomepred.tsv.gz" , paste(ugene_name2,":",pos_start,"-",pos_start,sep=""), " >> output/output.csv", sep =" ")
    #         }else if(input$endPosition=="" && input$startPosition==""){
    #           command_filtrationh <- paste("zcat data/domerefdata/final/proteomedomepred.tsv.gz | head -1 > output/output.csv")  
    #           command_filtration2 <- paste("tabix data/domerefdata/final/proteomedomepred.tsv.gz" , ugene_name2, " >> output/output.csv", sep =" ")
    #         } else{
    #           command_filtrationh <- paste("zcat data/domerefdata/final/proteomedomepred.tsv.gz | head -1 > output/output.csv")
    #           command_filtration2 <- paste("tabix data/domerefdata/final/proteomedomepred.tsv.gz" , paste(ugene_name2,":",pos_start,"-",pos_end,sep=""), " >> output/output.csv", sep =" ")
    #         }
    #         
    #         tier_sort = c(1,2,3,4,0)
    #         print(command_filtrationh)
    #         system(command_filtrationh, intern = TRUE)
    #         print(command_filtration2)
    #         system(command_filtration2, intern = TRUE)
    #         domain2 <- read.table(file = 'output/output.csv', header = T, sep = "\t")
    #         #print(domain2$tier)
    #         #domain2_sort <- domain2[match(tier_sort,domain2$tier),]
    #         domain2_1 <-subset(domain2,domain2$tier == 1)
    #         domain2_2 <-subset(domain2,domain2$tier == 2)
    #         domain2_3 <-subset(domain2,domain2$tier == 3)
    #         domain2_4 <-subset(domain2,domain2$tier == 4)
    #         domain2_0 <-subset(domain2,domain2$tier == 0)
    #         domain2_sort <- rbind(domain2_1,domain2_2,domain2_3,domain2_4,domain2_0)
    #         domain2_sort <- domain2_sort[ -c(1,14,15,16,18) ]
    #         #d <- mutate(domain2_sort, uniprot = ifelse((domain2_sort$"uniprotsite" != "-") & (domain2_sort$"X3d_uniprotsite_close" == "-"), domain2_sort$"uniprotsite",ifelse((domain2_sort$"uniprotsite" == "-") & (domain2_sort$"X3d_uniprotsite_close" != "-"), domain2_sort$"X3d_uniprotsite_close","-")))
    #         #d <- mutate(d, phosphosite_new = ifelse((d$"phosphosite" != "-") & (d$"X3d_pspsite_close" == "-"), d$"phosphosite",ifelse((d$"phosphosite" == "-") & (d$"X3d_pspsite_close" != "-"), d$"X3d_pspsite_close","-")))
    #         #d <- mutate(domain2_sort, uniprot = ifelse((domain2_sort$"uniprotsite" != "-") & (domain2_sort$"X3d_uniprotsite_close" == "-"), domain2_sort$"uniprotsite",ifelse((domain2_sort$"uniprotsite" == "-") & (domain2_sort$"X3d_uniprotsite_close" != "-"), domain2_sort$"X3d_uniprotsite_close",ifelse((domain2_sort$"phosphosite" != "-") & (domain2_sort$"X3d_pspsite_close" == "-"), domain2_sort$"phosphosite",ifelse((domain2_sort$"phosphosite" == "-") & (domain2_sort$"X3d_pspsite_close" != "-"), domain2_sort$"X3d_pspsite_close","-")))))
    #         #d <- mutate(d, hotspot = ifelse((d$"hotspot3dann" != "-") & (d$"X3d_hotspot_close" == "-"), d$"hotspot3dann",ifelse((d$"hotspot3dann" == "-") & (d$"X3d_hotspot_close" != "-"), d$"X3d_hotspot_close","-")))
    #         #d <- mutate(d, analog = ifelse((d$"X3d_statsig_close" != "-") & (d$"X3d_analog_close" == "-"), d$"X3d_statsig_close",ifelse((d$"X3d_statsig_close" == "-") & (d$"X3d_analog_close" != "-"), d$"X3d_analog_close","-")))
    #         #domain_final1 <- domain2_sort[-c(12,13,14,15,16,17,18,19)]
    #         #print(domain_final1)
    #         #colnames(domain_final1) <- c("Tier","Protein","Position","RefAA","AltAA","Domain","Hotspot/Resistant Mutation","Analogous to mutation","CADD score","COSMIC count","Clinvar","Interface","Functional site (Uniprot)","Functional site (Phosphosite)","Known 3D hotspot","Hotspot/analog proximity")
    #         #colnames(domain_final1) <- c("Tier","Protein","Position","RefAA","AltAA","Domain","Hotspot/Resistant Mutation","Analogous to mutation","CADD score","COSMIC count","Clinvar","Interface","Functional site","Known 3D hotspot","Hotspot/analog proximity")
    #         #print(domain2)
    #         output$domain_table2 <- renderDataTable(expr = domain2_sort, options = list(scrollX=T)
    #         )
    #         # Downloadable csv of selected dataset ----
    #         observe({
    #             shinyjs::toggleState(id = "downloadData2", !is.null(domain2_sort))
    #         })
    #         output$downloadData2 <- downloadHandler(
    #             filename = function() {
    #                 paste("output_rohit1", ".csv", sep = "")
    #             },
    #             content = function(file) {
    #                 write.csv(domain2, file, row.names = FALSE)
    #             }
    #         )
    #         Sys.sleep(2)
    #         
    #         
    #         #df <- as.data.frame(do.call(cbind, list(filname,variant_count)))
    #         #colnames(df) <- c("Input Files", "Varaints Counts")
    #         #print(df)
    #         #output$variants_count <- renderDataTable(expr = df)
    #         shinyjs::reset("form2")
    #         shinyjs::hide("form2")
    #         shinyjs::show("thankyou_msg2")
    #     },
    #     error = function(err) {
    #         shinyjs::html("error_msg2", err$message)
    #         shinyjs::show(id = "error2", anim = TRUE, animType = "fade")
    #     },
    #     finally = {
    #         shinyjs::enable("submit2")
    #         shinyjs::hide("submit_msg2")
    #     })
    # })
    
}

# Run the application
shinyApp(ui = ui, server = server)
