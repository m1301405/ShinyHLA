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

# If you need to ensure that when the second selection imgthla is updated, 
# the third selection package is also automatically updated, you can handle it this way
observeEvent(input$imgthla, {
  req(input$imgthla)
  updateSelectInput(session, "package", choices = filtered_package())
})

tmp_dir <- tempdir()  # Get the temporary directory

observeEvent(input$step2, {
  # Define file paths
  file1 <- file.path(tmp_dir, "Input", "sample1_1.fastq.gz")
  file2 <- file.path(tmp_dir, "Input", "sample1_2.fastq.gz")
  
  # Check if files exist
  if (file.exists(file1) && file.exists(file2)) {
    updateTabItems(session, "my_tabs", "hla")
  } else {
    # If files do not exist, display a simple warning message
    shinyalert("Error", "Please upload the sequencing file before the HLA typing step.", type = "error")
  }
})