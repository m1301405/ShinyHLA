## Fastq Upload & Warning
observeEvent(input$fastq_upload_button, {
  tmp_dir <- tempdir()  # Get the temporary directory
  bam_dir <- file.path(tmp_dir, "BAM")  # Define the BAM directory path
  input_dir <- file.path(tmp_dir, "Input")  # Define the Input directory path
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # Transparent background
  )
  w$show() # Show the loading animation
  
  if (is.null(input$fastq_gz_upload)) {
    output$file_count <- renderText("No files selected")
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
  }
  
  required_files <- 2
  num_files <- nrow(input$fastq_gz_upload)
  
  if (num_files != required_files) {
    output$file_count <- renderText(
      paste("You have selected", num_files, "files. Exactly", required_files, "files are required.")
    )
    shinyalert("Error", paste("You have selected", num_files, "files. Exactly", required_files, "files are required."), type = "error")
  } else {
    output$file_count <- renderText(
      paste("You have selected", num_files, "files. Upload successful.")
    )
    shinyalert("Success", paste("You have selected", num_files, "files. Upload successful."), type = "success")
    
    # Remove existing BAM directory and recreate it
    if (dir.exists(bam_dir)) {
      unlink(bam_dir, recursive = TRUE)
    }
    dir.create(bam_dir, showWarnings = FALSE)
    
    # Remove existing Input directory and recreate it
    if (dir.exists(input_dir)) {
      unlink(input_dir, recursive = TRUE)
    }
    dir.create(input_dir, showWarnings = FALSE)
    
    # Copy uploaded files to the Input directory
    file.copy(input$fastq_gz_upload$datapath, input_dir, overwrite = TRUE)
    
    # Reset the file input component
    session$sendCustomMessage(type = 'resetFileInput', message = list(id = "fastq_gz_upload"))
    
    # Rename uploaded files with meaningful names
    system(paste("mv", file.path(input_dir, "0.gz"), file.path(input_dir, "sample1_1.fastq.gz")))
    system(paste("mv", file.path(input_dir, "1.gz"), file.path(input_dir, "sample1_2.fastq.gz")))
  }
  w$hide()
})

## BAM Upload & Warning & Convert
observeEvent(input$bam_upload_button, {
  tmp_dir <- tempdir()  # Get the temporary directory
  bam_dir <- file.path(tmp_dir, "BAM")  # Define the BAM directory path
  input_dir <- file.path(tmp_dir, "Input")  # Define the Input directory path
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # Transparent background
  )
  w$show() # Show the loading animation
  
  if (is.null(input$bam_upload)) {
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
  }
  
  # Remove existing BAM directory and recreate it
  if (dir.exists(bam_dir)) {
    unlink(bam_dir, recursive = TRUE)
  }
  dir.create(bam_dir, showWarnings = FALSE)
  
  # Copy the uploaded BAM file to the BAM directory
  file.copy(input$bam_upload$datapath, bam_dir, overwrite = TRUE)
  
  # Convert BAM to Fastq
  bam_file <- file.path(bam_dir, "0.bam")
  system(paste("/arcasHLA/arcasHLA extract", bam_file, "-o", bam_dir, "-t", num_cores, "-v"))
  
  # Remove existing Input directory and recreate it
  if (dir.exists(input_dir)) {
    unlink(input_dir, recursive = TRUE)
  }
  dir.create(input_dir, showWarnings = FALSE)
  
  # Move sorted files to the Input directory
  system(paste("mv", file.path(bam_dir, "0.sorted.1.fq.gz"), file.path(input_dir, "sample1_1.fastq.gz")))
  system(paste("mv", file.path(bam_dir, "0.sorted.2.fq.gz"), file.path(input_dir, "sample1_2.fastq.gz")))
  
  shinyalert("Success", paste("Bam file Upload successful."), type = "success")
  
  w$hide()
})
