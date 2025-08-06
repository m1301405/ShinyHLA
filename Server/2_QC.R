## QC
# Custom function: Convert Coverage(%) into a pie chart with numeric values
coverage_to_pie <- function(value) {
  paste0(
    "<div style='
      position: relative;
      width: 50px;
      height: 50px;
      background: conic-gradient(
        orange ", value, "%,
        lightgray ", value, "%);
      border-radius: 50%;
      display: inline-block;'>
      <span style='
        position: absolute;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        font-size: 14px;
        font-weight: bold;
        color: black;'>
        ", value, "%
      </span>
    </div>"
  )
}
#====================================
observeEvent(input$coverage_button, {
  tmp_dir <- tempdir()
  qc_dir <- file.path(tmp_dir, "QC")
  #====================================
  bam_dir <- file.path(tmp_dir, "BAM")
  bam_file <- list.files(bam_dir, full.names = TRUE)
  
  if (length(bam_file) == 0) {
    shinyalert(
      title = "Error",
      text = "No input files found. Please upload file before QC",
      type = "error"
    )
    return()
  }
  #====================================
  # Create waiting animation
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # transparent background
  )
  w$show() # show the waiter

  #====================================
  # Calculate RNA-seq coverage
  ## Define variables for input file, gene annotation file, and output directory
  
  input_bam <- file.path(tmp_dir, "BAM", "sample.bam")
  annotation_gff <- "/root/shiny/PanDepth/hg38/gencode.v38.annotation.gff3"
  output_dir <- file.path(tmp_dir, "QC", "hla_coverage")
  

  ## QC folder
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## Construct command
  command <- sprintf(
    "/root/shiny/PanDepth/bin/pandepth -i %s -g %s -o %s",
    input_bam, annotation_gff, output_dir
  )

  ## Execute command
  system(command)

  #====================================
  # Decompress
  ## Define the compressed file and the target decompression path
  
  input_gz <- file.path(tmp_dir, "QC", "hla_coverage.gene.stat.gz")
  output_file <- file.path(tmp_dir, "QC", "hla_coverage.gene.stat")

  if (file.exists(output_file)) {
    file.remove(output_file)
  }

  ## Decompress file
  gunzip(input_gz, destname = output_file, remove = FALSE)
  #====================================
  # Extract HLA region data
  ## Define the path of the Shell script
  qc_pipeline_script <- "/root/shiny/Server/hla_qc_pipline.sh"
  output_file_coverage <- file.path(tmp_dir, "QC", "hla_coverage.csv")
  ## Execute Shell script, passing output_file as a parameter
  command <- sprintf("bash %s %s %s", qc_pipeline_script, output_file, output_file_coverage)
  system(command)

  #====================================
  # Generate coverage output table
  output$coverage_table <- renderDataTable({
    # Load the output CSV file
    hla_coverage <- read_delim(
      output_file_coverage,
      delim = "\t", escape_double = FALSE, trim_ws = TRUE
    )

    # Format the table
    formatted_table <- hla_coverage %>%
      mutate(
        `Coverage(%)` = sapply(`Coverage(%)`, coverage_to_pie)
      )

    # Generate DataTable
    datatable(
      formatted_table,
      options = list(
        autoWidth = TRUE,
        scrollX = TRUE,
        responsive = TRUE,
        deferRender = TRUE, # Delay rendering to prevent misalignment
        pageLength = 10,
        columnDefs = list(
          list(width = '10%', targets = 0), # Explicitly specify the width of each column
          list(width = '15%', targets = 1:2),
          list(width = '20%', targets = 3:4),
          list(className = "dt-center", targets = "_all") # Alignment
        ),
        drawCallback = JS(
          "function(settings) {",
          "  this.api().columns.adjust();", # Force column width adjustment
          "  $(window).trigger('resize');", # Trigger window resize event
          "}"
        )
      ),
      rownames = FALSE, # Remove rownames
      escape = FALSE,   # Support HTML format
      class = "display nowrap stripe"
    )
  })

  # Hide waiting animation
  w$hide()

  # Add CSV download functionality
  output$coverage_download <- downloadHandler(
    filename = function() {
      paste0("HLA_Coverage_", Sys.Date(), ".csv") # Retain CSV extension
    },
    content = function(file) {
      # Read the original file
      data <- read_delim(
        output_file_coverage,
        delim = "\t", # Specify tab-delimited
        escape_double = FALSE,
        trim_ws = TRUE
      )
      
      # Write data as comma-separated CSV
      write_csv(data, file)
    }
  )
})

