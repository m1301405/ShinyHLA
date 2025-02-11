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
                           choices = list("FASTQ" = "fastq", "BAM" = "bam","Demo" = "demo"), inline = TRUE),
              conditionalPanel(
                condition = "input.upload_choice == 'fastq'",
                fileInput(inputId = "fastq_gz_upload", 
                          label = "",
                          accept = ".gz", 
                          multiple = TRUE, 
                          placeholder = "Please upload 2 paired-end FASTQ (.fastq.gz) files, with a maximum combined size of 5GB."),
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
                          placeholder = "Upload BAM file; maximum size: 5GB."),
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
                      title = "The demo uses HG00160 from 1000 Genomes Project (GRCh38), with chr6 reads extracted.",
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
              conditionalPanel(
                condition = "((input.bam_upload_button && input.upload_choice === 'bam') || (input.demo_upload_button && input.upload_choice === 'demo'))",
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
                      title = "Calculate Total Depth, Coverage, and Mean Depth for each HLA allele",
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
                style = "display: flex; justify-content: space-between; align-items: center; padding: 10px;",
                actionButton(
                  inputId = "step2",
                  label = "Next step: HLA typing",
                  icon = icon("arrow-circle-right"))
              )
            ),
            tags$div(style = "height: 20px;"),
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