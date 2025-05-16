## Fastq Upload & Warning
observeEvent(input$fastq_upload_button, {
  tmp_dir <- tempdir()
  bam_dir <- file.path(tmp_dir, "BAM")
  input_dir <- file.path(tmp_dir, "Input")
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  if (is.null(input$fastq_gz_upload)) {
    output$file_count <- renderText("No files selected")
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
  }
  
  required_files <- 2
  num_files <- nrow(input$fastq_gz_upload)
  
  if (num_files != required_files) {
    shinyalert("Error", paste("Please select exactly", required_files, "files. You selected", num_files, "."), type = "error")  } else {
  } else {
    if (dir.exists(bam_dir)) {
      unlink(bam_dir, recursive = TRUE)
    }
    dir.create(bam_dir, showWarnings = FALSE)
    
    if (dir.exists(input_dir)) {
      unlink(input_dir, recursive = TRUE)
    }
    dir.create(input_dir, showWarnings = FALSE)
    
    # Copy uploaded files to the Input directory
    file.copy(input$fastq_gz_upload$datapath, input_dir, overwrite = TRUE)
    
    # Rename uploaded files with meaningful names
    system(paste("mv", file.path(input_dir, "0.gz"), file.path(input_dir, "sample_1.fastq.gz")))
    system(paste("mv", file.path(input_dir, "1.gz"), file.path(input_dir, "sample_2.fastq.gz")))
    
    shinyalert("Success", paste("FASTQ Files Upload Successful."), type = "success")
  }
  w$hide()
})

## BAM Upload & Warning & Convert
observeEvent(input$bam_upload_button, {
  tmp_dir <- tempdir()
  bam_dir <- file.path(tmp_dir, "BAM")
  input_dir <- file.path(tmp_dir, "Input")
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # Transparent background
  )
  w$show() # Show the loading animation
  
  if (is.null(input$bam_upload)) {
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
    
  } else {
    if (dir.exists(bam_dir)) {
      unlink(bam_dir, recursive = TRUE)
    }
    dir.create(bam_dir, showWarnings = FALSE)
    
    file.copy(input$bam_upload$datapath, bam_dir, overwrite = TRUE)
    
    # Convert BAM to Fastq
    bam_file <- file.path(bam_dir, "sample.bam")
    system(paste("/arcasHLA/arcasHLA extract", bam_file, "-o", bam_dir, "-t", num_cores, "-v"))
    
    # Remove existing Input directory and recreate it
    if (dir.exists(input_dir)) {
      unlink(input_dir, recursive = TRUE)
    }
    dir.create(input_dir, showWarnings = FALSE)
    
    # Move sorted files to the Input directory
    system(paste("mv", file.path(bam_dir, "sample.sorted.1.fq.gz"), file.path(input_dir, "sample_1.fastq.gz")))
    system(paste("mv", file.path(bam_dir, "sample.sorted.2.fq.gz"), file.path(input_dir, "sample_2.fastq.gz")))
    
    shinyalert("Success", paste("Bam File Upload Successful."), type = "success")
  }
  w$hide()
})

## Demo QC option
wes_demo_uploaded_once <- reactiveVal(FALSE)

is_current_demo_wes <- reactive({
  input$upload_choice == "demo" && input$demo_choice == "WES"
})

observeEvent(input$demo_upload_button, {
  if (input$upload_choice == "demo") {
    if (input$demo_choice == "WES") {
      wes_demo_uploaded_once(TRUE)  
    } else if (input$demo_choice == "RNA-seq") {
      wes_demo_uploaded_once(FALSE)
    }
  }
})

output$show_qc_button <- reactive({
  (input$upload_choice == "bam" && input$bam_upload_button >= 1) ||
    (input$upload_choice == "demo" && is_current_demo_wes() && wes_demo_uploaded_once())
})
outputOptions(output, "show_qc_button", suspendWhenHidden = FALSE)

