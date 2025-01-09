filtered_imagthla <- reactive({
  req(input$sequence)
  HLA_typing %>%
    filter(NGS_type %in% input$sequence) %>%
    pull(IMGTHLA_version) %>%
    unique() %>%
    sort(decreasing = TRUE)
})

filtered_package <- reactive({
  req(input$sequence, input$imgthla)
  HLA_typing %>%
    filter(NGS_type %in% input$sequence, IMGTHLA_version %in% input$imgthla) %>%
    pull(Typing_package) %>%
    unique()
})

observe({
  req(input$sequence)
  updateSelectInput(session, "imgthla", choices = filtered_imagthla())
})

observe({
  req(input$sequence, input$imgthla)
  updateSelectInput(session, "package", choices = filtered_package())
})

# 如果需要確保第二個選項 imgthla 更新時，第三個選項 package 也會自動更新，可以這樣處理
observeEvent(input$imgthla, {
  req(input$imgthla)
  updateSelectInput(session, "package", choices = filtered_package())
})

tmp_dir <- tempdir()  # 获取临时目录

observeEvent(input$step2, {
  # 定义文件路径
  file1 <- file.path(tmp_dir, "Input", "sample1_1.fastq.gz")
  file2 <- file.path(tmp_dir, "Input", "sample1_2.fastq.gz")
  
  # 检查文件是否存在
  if (file.exists(file1) && file.exists(file2)) {
    updateTabItems(session, "my_tabs", "hla")
  } else {
    # 如果文件不存在，显示简单的警告信息
    shinyalert("Error", "Please upload the sequencing file before the HLA typing step.", type = "error")
  }
})

