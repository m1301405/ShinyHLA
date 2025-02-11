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
    # 移除現有BAM及其內容
    if (dir.exists("./BAM")) {
      unlink("./BAM", recursive = TRUE)
    }
    # 創建BAM並加入文件
    dir.create("./BAM", showWarnings = FALSE)
    # 複製BAM file
    file.copy("/root/shiny/Demo/wes/HG00160_wes.bam", "./BAM/0.bam", overwrite = TRUE)
    
    # 移除現有Input及其內容
    if (dir.exists("./Input")) {
      unlink("./Input", recursive = TRUE)
    }
    # 創建Input並覆蓋舊文件
    dir.create("./Input", showWarnings = FALSE)
    # copy Fastq
    file.copy("/root/shiny/Demo/wes/HG00160_wes_1.fq.gz", "./Input/sample1_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/wes/HG00160_wes_2.fq.gz", "./Input/sample1_2.fastq.gz", overwrite = TRUE)
    
    # alert
    shinyalert("Success", paste("WES Demo File Upload successful."), type = "success")
    
  } else if (input$demo_choice == "RNA-seq")  {
    # 移除現有BAM及其內容
    if (dir.exists("./BAM")) {
      unlink("./BAM", recursive = TRUE)
    }
    # 創建BAM並加入文件
    dir.create("./BAM", showWarnings = FALSE)
    # 複製BAM file
    file.copy("/root/shiny/Demo/rna-seq/HG00160_rseq.bam", "./BAM/0.bam", overwrite = TRUE)
    
    # 移除現有Input及其內容
    if (dir.exists("./Input")) {
      unlink("./Input", recursive = TRUE)
    }
    # 創建Input並覆蓋舊文件
    dir.create("./Input", showWarnings = FALSE)
    # copy Fastq
    file.copy("/root/shiny/Demo/rna-seq/HG00160_rseq_1.fq.gz", "./Input/sample1_1.fastq.gz", overwrite = TRUE)
    file.copy("/root/shiny/Demo/rna-seq/HG00160_rseq_2.fq.gz", "./Input/sample1_2.fastq.gz", overwrite = TRUE)
    
    # alert
    shinyalert("Success", paste("RNA-seq Demo File Upload successful."), type = "success")
  }
  w$hide() 
  
  setwd(ori_dir)
})

