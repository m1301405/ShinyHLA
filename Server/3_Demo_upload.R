observeEvent(input$demo_upload_button, {
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  if (input$demo_choice == "WES") {
    if (dir.exists("./BAM")) unlink("./BAM", recursive = TRUE)
    dir.create("./BAM", showWarnings = FALSE)
    file.copy("/root/shiny/Demo/wes/HG00096.bam", "./BAM/sample.bam", overwrite = TRUE)
    
    if (dir.exists("./Input")) unlink("./Input", recursive = TRUE)
    dir.create("./Input", showWarnings = FALSE)
    file.copy("/root/shiny/Demo/wes/HG00096_1.fastq.gz", "./Input/sample_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/wes/HG00096_2.fastq.gz", "./Input/sample_2.fastq.gz", overwrite = TRUE)
    
    source_mode("demo")
    demo_uploaded(TRUE)
    demo_selected_ngs(input$demo_choice)
    
    updateSelectInput(
      session,
      "sequence",
      choices = c(input$demo_choice),
      selected = input$demo_choice
    )
    
    shinyalert("Success", "WES Demo File Upload Successful.", type = "success")
    
  } else if (input$demo_choice == "RNA-seq") {
    if (dir.exists("./BAM")) unlink("./BAM", recursive = TRUE)
    dir.create("./BAM", showWarnings = FALSE)
    
    if (dir.exists("./Input")) unlink("./Input", recursive = TRUE)
    dir.create("./Input", showWarnings = FALSE)
    file.copy("/root/shiny/Demo/rna-seq/HG00096_1.fastq.gz", "./Input/sample_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/rna-seq/HG00096_2.fastq.gz", "./Input/sample_2.fastq.gz", overwrite = TRUE)
    
    source_mode("demo")
    demo_uploaded(TRUE)
    demo_selected_ngs(input$demo_choice)
    
    updateSelectInput(
      session,
      "sequence",
      choices = c(input$demo_choice),
      selected = input$demo_choice
    )
    
    shinyalert("Success", "RNA-seq Demo File Upload Successful.", type = "success")
  }
  
  w$hide()
  setwd(ori_dir)
})

