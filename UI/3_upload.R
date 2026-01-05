tabItem(tabName = "upload",
        fluidRow(
          column(
            width = 12,
            box(
              id = "file_upload_box",
              height = "auto",
              width = 12,
              elevation = 4,
              title = strong("Data Upload"),
              solidHeader = TRUE,
              maximizable = TRUE,
              status = "primary",
              radioButtons("upload_choice", "Choose data format:",
                           choices = list("FASTQ" = "fastq", "BAM" = "bam","Demo" = "demo", "JobID" = "jobid"), inline = TRUE),
              conditionalPanel(
                condition = "input.upload_choice == 'fastq'",
                fileInput(inputId = "fastq_gz_upload", 
                          label = "",
                          accept = ".gz", 
                          multiple = TRUE, 
                          placeholder = "Upload 2 paired-end fastq.gz files; Total max: 5GB"),
                actionButton(
                  inputId = "fastq_upload_button", 
                  label = "Upload",
                  icon = icon("upload")
                )
              ),
              conditionalPanel(
                condition = "input.upload_choice == 'bam'",
                fileInput(inputId = "bam_upload", 
                          label = "",
                          accept = ".bam", 
                          placeholder = "Upload bam file; max: 5GB"),
                actionButton(
                  inputId = "bam_upload_button", 
                  label = "Upload",
                  icon = icon("upload")
                )
              ),
              conditionalPanel(
                condition = "input.upload_choice == 'demo'",
                awesomeRadio(
                  inputId = "demo_choice",
                  label = tagList(
                    "Demo File", 
                    tooltip(
                      icon("info-circle"),
                      title = "Demo uses HG00096 data: WES (aligned BAM) and RNA-Seq (raw FASTQ) for HLA typing.",
                      placement = "right"
                    )
                  ),
                  choices = c("WES", "RNA-seq"),
                  selected = "WES",
                  checkbox = TRUE
                ),
                actionButton(
                  inputId = "demo_upload_button", 
                  label = "Upload",
                  icon = icon("upload")
                )
              ),
              # 2. [新增] JobID Panel
              conditionalPanel(
                condition = "input.upload_choice == 'jobid'",
                div(
                  style = "display: flex; align-items: center; margin-bottom: 10px;",
                  div(
                    style = "flex-grow: 1; margin-right: 10px;",
                    textInput(
                      inputId = "jobid_input", 
                      label = NULL, 
                      placeholder = "Enter JobID (e.g., 20260105-ABCD)"
                    )
                  ),
                  actionButton(
                    inputId = "jobid_load_button", 
                    label = "Load Result",
                    icon = icon("sync-alt"),
                    class = "btn-info" # 使用不同顏色區分功能
                  )
                ),
                tags$p("Load previously analyzed results directly without re-running the analysis.", 
                       style = "font-size: 12px; color: #666;")
              ),
              conditionalPanel(
                condition = "output.show_qc_button",
                tags$hr(style = "border-top: 2px solid #000;"),
                div(
                  style = "display: flex; align-items: center;",
                  actionButton(
                    inputId = "coverage_button",
                    label = "Quality control",
                    icon = icon("chart-line")
                  ),
                  div(
                    style = "margin-left: 10px;",
                    tooltip(
                      icon("info-circle"),
                      title = "Check total depth, coverage, and mean depth per HLA allele to ensure data quality.",
                      placement = "right"
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
            conditionalPanel(
              condition = "input.fastq_upload_button >= 1 || input.bam_upload_button >= 1 || input.demo_upload_button >=1",
              div(
                style = "display: flex; align-items: center; padding: 10px;",
                div(
                  style = "display: flex; align-items: center;",
                  actionButton(
                    inputId = "step2",
                    label = "Next step: HLA typing",
                    icon = icon("arrow-circle-right")
                  ),
                  div(
                    style = "margin-left: 8px;",
                    tooltip(
                      icon("exclamation-triangle"),
                      title = "Only the most recently uploaded file is retained for analysis. To re-analyze a previously uploaded file, please re-upload it.",
                      placement = "right")
                    )
                  )
                )
              ),
            tags$div(style = "height: 20px;")
            )
          ),
        fluidRow(
          column(
            width = 12,
            conditionalPanel(
              condition = "input.coverage_button >= 1",
              box(
                id = "QC_box",
                height = "auto",
                width = 12,
                elevation = 4,
                title = strong("Sequencing coverage"),
                solidHeader = FALSE,
                maximizable = TRUE,
                status = "info",
                dataTableOutput("coverage_table")
              )
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            conditionalPanel(
              condition = "input.coverage_button >= 1",
              div(
                style = "display: flex; justify-content: space-between; align-items: center; padding: 10px;",
                downloadButton(
                  outputId = "coverage_download",
                  label = "Download")
              )
            )
          )
        )
)