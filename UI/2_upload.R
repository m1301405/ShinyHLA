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
                          placeholder = "Upload 2 paired-end fastq.gz files; max:15GB"),
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
                          placeholder = "Upload bam file; max:15GB"),
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
                      title = "The demo uses HG00096 from the 1000 Genomes project on the GRCh38 dataset, pre-processed by extracting chr6 reads.",
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
                      title = "Calculate Total Depth, Coverage, and Mean Depth for each HLA allele to ensure high-quality input data and reliable results.",
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