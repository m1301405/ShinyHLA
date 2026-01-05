observeEvent(input$demo_upload_button, {
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # transparent background
  )
  w$show() # show the waiter
  
  # Demo upload 
  if (input$demo_choice == "WES") {
    # Remove existing BAM and its contents
    if (dir.exists("./BAM")) {
      unlink("./BAM", recursive = TRUE)
    }
    # Create BAM directory and add files
    dir.create("./BAM", showWarnings = FALSE)
    # Copy BAM file
    file.copy("/root/shiny/Demo/wes/HG00096.bam", "./BAM/sample.bam", overwrite = TRUE)
    
    # Remove existing Input and its contents
    if (dir.exists("./Input")) {
      unlink("./Input", recursive = TRUE)
    }
    # Create Input directory and overwrite old files
    dir.create("./Input", showWarnings = FALSE)
    # Copy Fastq files
    file.copy("/root/shiny/Demo/wes/HG00096_1.fastq.gz", "./Input/sample_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/wes/HG00096_2.fastq.gz", "./Input/sample_2.fastq.gz", overwrite = TRUE)
    
    # Alert
    shinyalert("Success", paste("WES Demo File Upload Successful."), type = "success")
    
  } else if (input$demo_choice == "RNA-seq")  {
    # Remove existing BAM and its contents
    if (dir.exists("./BAM")) {
      unlink("./BAM", recursive = TRUE)
    }
    # Create BAM directory and add files
    dir.create("./BAM", showWarnings = FALSE)

    # Remove existing Input and its contents
    if (dir.exists("./Input")) {
      unlink("./Input", recursive = TRUE)
    }
    # Create Input directory and overwrite old files
    dir.create("./Input", showWarnings = FALSE)
    # Copy Fastq files
    file.copy("/root/shiny/Demo/rna-seq/HG00096_1.fastq.gz", "./Input/sample_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/rna-seq/HG00096_2.fastq.gz", "./Input/sample_2.fastq.gz", overwrite = TRUE)
    
    # Alert
    shinyalert("Success", paste("RNA-seq Demo File Upload Successful."), type = "success")
  }
  w$hide() 
  
  setwd(ori_dir)
})