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
                  selectInput("sequence",label = "NGS type:", choices = HLA_typing %>% pull(NGS_type) %>% unique(), selected = "WES"),
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
                      tooltip(
                        icon("info-circle"),
                        title = "Select the sequencing type, HLA reference genome version, and genotyping tool for analysis.",
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
                      tooltip(
                        icon("info-circle"),
                        title = "Use the MHC class I & II buttons to organize genotyping results into tables. Adjust columns and rows by dragging and dropping items.",
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
                      fluidRow(
                        conditionalPanel(
                          condition = "input.hla_typing_button >= 1",
                          uiOutput("summary_boxes"),
                          DTOutput("hla_typing_table",width = "auto",height = "auto"),
                          tags$br()
                        ),
                        div(
                          style = "display: flex; justify-content: space-between; align-items: center; width: 100%;",
                          conditionalPanel(
                            condition = "input.hla_typing_button >= 1 && input.package_value_ui != 'arcasHLA'",
                            div(
                              style = "flex-grow: 1;",
                              actionButton(
                                inputId = "igv_reference",
                                label = "IGV",
                                icon = icon("dna")
                              ),
                              tooltip(
                                icon("info-circle"),
                                title = "Use the Integrative Genomics Viewer (IGV) to visualize the alignment results from the typing process.",
                                placement = "right"
                              )
                            )
                          ),
                          conditionalPanel(
                            condition = "input.hla_typing_button >= 1",
                            div(
                              style = "flex-grow: 0;",
                              downloadButton(outputId = "hla_typing_download", label = "Download"),
                              tooltip(
                                icon("info-circle"),
                                title = "Save the nucleotide and protein sequences of HLA alleles in fasta format.",
                                placement = "right"
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
                                labels = list("NGS_type" = "NGS Type", 
                                              "Tool" = "Tool"),
                                input_id = "mhci_col_list"
                              ),
                              add_rank_list(
                                text = HTML("Rows (<i>maximum=2</i>)"),
                                labels = list("MHC_class" = "MHC Class",
                                              "IMGTHLA_version" = "IMGTHLA Version"),
                                input_id = "mhci_row_list"
                              )
                            )
                          )
                        ),
                        column(
                          width = 12,
                          uiOutput("mhci_warning_message")
                        ),
                        column(
                          width = 12,
                          div(class = "custom-box", style = "overflow-x: auto;", pivottablerOutput("mhci_pivot_table")
                          )
                        ),
                        column(
                          width = 12,
                          tags$br(),
                          conditionalPanel(
                            condition = "input.mhci_pivottable_button >= 1",
                            downloadButton(
                              outputId = "mhci_pivot_download",
                              label = "Download"
                            )
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
                                labels = list("NGS_type" = "NGS Type", 
                                              "Tool" = "Tool"),
                                input_id = "mhcii_col_list"
                              ),
                              add_rank_list(
                                text = HTML("Rows (<i>maximum=2</i>)"),
                                labels = list("MHC_class" = "MHC Class",
                                              "IMGTHLA_version" = "IMGTHLA Version"),
                                input_id = "mhcii_row_list"
                              )
                            )
                          )
                        ),
                        column(
                          width = 12,
                          uiOutput("mhcii_warning_message")
                        ),
                        column(
                          width = 12,
                          div(class = "custom-box", style = "overflow-x: auto;", pivottablerOutput("mhcii_pivot_table")
                          )
                        ),
                        column(
                          width = 12,
                          tags$br(),
                          conditionalPanel(
                            condition = "input.mhcii_pivottable_button >= 1",
                            downloadButton(
                              outputId = "mhcii_pivot_download",
                              label = "Download")
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
          column(
            width = 12,
            box(
              id = "igv_box",
              height = "auto",
              width = 12,
              elevation = 4,
              title = strong("Integrative Genomics Viewer"),
              collapsible = TRUE,
              solidHeader = FALSE,
              maximizable = TRUE,
              status = "danger",
              hidden(
                textInput("package_value_ui", label = NULL, value = "")
              ),
              div(
                style = "border-left: 4px solid #d9534f; padding: 10px; background-color: #f9f9f9; margin-bottom: 15px;",
                HTML(
                  "<p style='font-size: 12px; color: #555; margin: 0;'>
           <strong>Hint:</strong> Due to the incomplete integration between the IGV and Shiny packages, 
           it is recommended to adjust the webpage window to an appropriate size before proceeding with 'Display'. 
           If you notice any anomalies in the Genome Track display, 
           you can correct the issue by readjusting the window size of the ShinyHLA webpage.
           </p>"
                )
              ),
              fluidRow(
                column(
                  width = 12,
                  tabsetPanel(
                    id = "igv_tabsetPanel",
                    tabPanel(
                      title = "GRCh38-I",
                      fluidRow(
                        column(
                          width = 12,
                          igvShinyOutput('hg38_igv_i', width = "100%", height = "600px")
                        )
                      ),
                      fluidRow(
                        column(
                          width = 12,
                          tags$br(),
                          div(
                            style = "display: flex; align-items: center;",
                            actionButton(
                              inputId = "hg38_igv_i_bam_button",
                              label = "Display",
                              icon = icon("eye"),
                            ),
                            div(
                              style = "margin-left: 10px;",
                              tooltip(
                                icon("info-circle"),
                                title = "If you provide sequencing files in BAM format, you can upload them for sequence visualization (GRCh38-aligned files only).",
                                placement = "right"
                              )
                            )
                          )
                        )
                      )
                    ),
                    tabPanel(
                      title = "GRCh38-II",
                      fluidRow(
                        column(
                          width = 12,
                          igvShinyOutput('hg38_igv_ii', width = "100%", height = "600px")
                        )
                      ),
                      fluidRow(
                        column(
                          width = 12,
                          tags$br(),
                          div(
                            style = "display: flex; align-items: center;",
                            actionButton(
                              inputId = "hg38_igv_ii_bam_button",
                              label = "Display",
                              icon = icon("eye"),
                            ),
                            div(
                              style = "margin-left: 10px;",
                              tooltip(
                                icon("info-circle"),
                                title = "If you provide sequencing files in BAM format, you can upload them for sequence visualization (GRCh38-aligned files only).",
                                placement = "right"
                              )
                            )
                          )
                        )
                      )
                    ),
                    tabPanel(
                      title = "Package",
                      conditionalPanel(
                        condition = "input.igv_reference >= 1",
                        fluidRow(
                          column(
                            width = 12,
                            igvShinyOutput("alignment_igv", width = "100%", height = "600px")
                          )
                        )
                      ),
                      conditionalPanel(
                        condition = "input.package_value_ui == 'HLA-HD'",
                        tags$br(),
                        fluidRow(
                          column(
                            width = 12,
                            igvShinyOutput("alignment_ii_igv", width = "100%", height = "600px")
                          )
                        )
                      ),
                      conditionalPanel(
                        condition = "input.igv_reference >= 1 && input.package_value_ui != 'arcasHLA'" ,
                        tags$br(),
                        fluidRow(
                          column(
                            width = 12,
                            div(
                              style = "display: flex; align-items: center;",
                              actionButton(
                                inputId = "igv_bam_button",
                                label = "Display Alignment File",
                                icon = icon("eye"),
                                style = "margin-right: 10px;"  
                              ),
                              div(
                                style = "margin-left: 10px;",
                                tooltip(
                                  icon("info-circle"),
                                  title = "Upload the BAM file generated during the HLA typing process after sequence alignment.",
                                  placement = "right"
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
)
