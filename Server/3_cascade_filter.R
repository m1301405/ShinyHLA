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

observeEvent(input$imgthla, {
  req(input$imgthla)
  updateSelectInput(session, "package", choices = filtered_package())
})

tmp_dir <- tempdir()  # Get the temporary directory

observeEvent(input$step2, {
  file1 <- file.path(tmp_dir, "Input", "sample_1.fastq.gz")
  file2 <- file.path(tmp_dir, "Input", "sample_2.fastq.gz")
  
  if (file.exists(file1) && file.exists(file2)) {
    updateTabItems(session, "my_tabs", "hla")
  } else {
    shinyalert(
      title = "Error",
      text = "Please upload the correct number and format of sequencing files before the HLA typing step.",
      type = "error"
    )
  }
})