tabItem(tabName = "hla",
        fluidRow(
          column(
            width = 12,
            box(
              id = "hla_box",
              height = "auto",
              width = 12,
              elevation = 4,
              title = strong("HLA Typing"),
              collapsible = TRUE,
              solidHeader = FALSE,
              maximizable = TRUE,
              status = "success",
              
              fluidRow(
                column(
                  width = 3,
                  
                  selectInput("sequence",label = "NGS type:", 
                              choices = HLA_typing %>% pull(NGS_type) %>% unique(), 
                              selected = "WES"),
                  
                  selectInput("imgthla", label = "IPD-IMGT/HLA version:", choices = NULL),
                  selectInput("package", label = "Typing tool:", choices = NULL),
                  
                  div(
                    style = "display: flex; align-items: center;",
                    actionButton(
                      inputId = "hla_typing_button", 
                      label = "Submit", 
                      icon = icon("rocket")
                    ),
                    div(
                      style = "margin-left: 10px;",
                      bs4Dash::tooltip(
                        icon("info-circle"),
                        title = "Select sequencing type, HLA reference genome version, and genotyping tool for analysis.",
                        placement = "right"
                      )
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "input.hla_typing_button >= 1",
                    
                    tags$hr(),
                    tags$p("Summarize:", style = "font-weight: bold; margin-right: 10px;"),
                    
                    div(
                      style = "display: flex; align-items: center;",
                      
                      actionButton(
                        inputId = "mhci_pivottable_button",
                        label = "MHC class I",
                        style = "margin-right: 10px;"
                      ),
                      
                      actionButton(
                        inputId = "mhcii_pivottable_button",
                        label = "MHC class II",
                        style = "margin-right: 10px;"
                      ),
                      
                      bs4Dash::tooltip(
                        icon("info-circle"),
                        title = "Use the MHC Class I & II buttons to organize genotyping results into tables.",
                        placement = "right"
                      )
                    ),
                    
                    tags$hr(),
                    
                    div(
                      style = "height: 200px; overflow-y: auto;",
                      uiOutput("hla_log")
                    )
                  )
                ),
                
                column(
                  width = 9,
                  
                  tabsetPanel(
                    id = "hla_panel",
                    
                    tabPanel(
                      title = "Result",
                      tags$br(),
                      
                      tags$div(
                        class = "hla-result-scroll",
                        tags$div(
                          class = "hla-result-inner",
                          
                          fluidRow(
                            conditionalPanel(
                              condition = "input.hla_typing_button >= 1 || input.jobid_load_button >= 1",
                              
                              uiOutput("summary_boxes"),
                              tags$hr(style = "margin: 20px 0;"),
                              
                              DTOutput("hla_typing_table", width = "100%"),
                              tags$br()
                            ),
                            
                            tags$div(
                              style = "display: flex; justify-content: flex-end; align-items: center; width: 100%; margin-top: 10px;",
                              
                              conditionalPanel(
                                condition = "input.hla_typing_button >= 1 || input.jobid_load_button >= 1",
                                div(
                                  downloadButton(
                                    outputId = "hla_typing_download", 
                                    label = "Download", 
                                    class = "btn-outline-success"
                                  ),
                                  bs4Dash::tooltip(icon("info-circle"), title = "Save sequences in FASTA format.", placement = "right")
                                )
                              )
                            )
                          )
                        )
                      )
                    ),
                    
                    tabPanel(
                      title = "MHC class I",
                      fluidRow(
                        column(
                          width = 12,
                          conditionalPanel(
                            condition = "input.mhci_pivottable_button >= 1",
                            bucket_list(
                              header = "",
                              group_name = "mhci_list",
                              orientation = "horizontal",
                              add_rank_list(
                                text = "Placement of unneeded items",
                                labels = NULL,
                                input_id = "mhci_ColRow_list"
                              ),
                              add_rank_list(
                                text = HTML("Columns (<i>maximum=2</i>)"),
                                labels = list("NGS_type" = "NGS Type","Tool" = "Tool"),
                                input_id = "mhci_col_list"
                              ),
                              add_rank_list(
                                text = HTML("Rows (<i>maximum=2</i>)"),
                                labels = list("MHC_class" = "MHC Class","IMGTHLA_version" = "IMGTHLA Version"),
                                input_id = "mhci_row_list"
                              )
                            )
                          )
                        ),
                        column(width = 12, uiOutput("mhci_warning_message")),
                        column(width = 12,
                               div(class = "custom-box", style = "overflow-x: auto;",
                                   pivottablerOutput("mhci_pivot_table"))
                        ),
                        column(width = 12,
                               tags$br(),
                               conditionalPanel(
                                 condition = "input.mhci_pivottable_button >= 1",
                                 downloadButton("mhci_pivot_download","Download")
                               )
                        )
                      )
                    ),
                    
                    tabPanel(
                      title = "MHC class II",
                      fluidRow(
                        column(
                          width = 12,
                          conditionalPanel(
                            condition = "input.mhcii_pivottable_button >= 1",
                            bucket_list(
                              header = "",
                              group_name = "mhcii_list",
                              orientation = "horizontal",
                              add_rank_list(
                                text = "Placement of unneeded items",
                                labels = NULL,
                                input_id = "mhcii_ColRow_list"
                              ),
                              add_rank_list(
                                text = HTML("Columns (<i>maximum=2</i>)"),
                                labels = list("NGS_type" = "NGS Type","Tool" = "Tool"),
                                input_id = "mhcii_col_list"
                              ),
                              add_rank_list(
                                text = HTML("Rows (<i>maximum=2</i>)"),
                                labels = list("MHC_class" = "MHC Class","IMGTHLA_version" = "IMGTHLA Version"),
                                input_id = "mhcii_row_list"
                              )
                            )
                          )
                        ),
                        column(width = 12, uiOutput("mhcii_warning_message")),
                        column(width = 12,
                               div(class = "custom-box", style = "overflow-x: auto;",
                                   pivottablerOutput("mhcii_pivot_table"))
                        ),
                        column(width = 12,
                               tags$br(),
                               conditionalPanel(
                                 condition = "input.mhcii_pivottable_button >= 1",
                                 downloadButton("mhcii_pivot_download","Download")
                               )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        ),
        fluidRow(
          div(
            id = "pharmgkb_container",
            uiOutput("pharmgkb_box_ui")
          )
        )
)