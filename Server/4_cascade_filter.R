# ============================================================
# Mode-aware cascade choices
# ============================================================
filtered_imagthla <- reactive({
  req(input$sequence)
  
  if (identical(source_mode(), "demo")) {
    demo_df <- read_demo_list()
    
    selected_ngs <- demo_selected_ngs()
    req(selected_ngs)
    
    demo_df %>%
      dplyr::filter(NGS_type %in% selected_ngs) %>%
      dplyr::pull(IMGTHLA_version) %>%
      unique() %>%
      sort(decreasing = TRUE)
  } else {
    HLA_typing %>%
      dplyr::filter(NGS_type %in% input$sequence) %>%
      dplyr::pull(IMGTHLA_version) %>%
      unique() %>%
      sort(decreasing = TRUE)
  }
})

filtered_package <- reactive({
  req(input$sequence, input$imgthla)
  
  if (identical(source_mode(), "demo")) {
    demo_df <- read_demo_list()
    
    selected_ngs <- demo_selected_ngs()
    req(selected_ngs)
    
    demo_df %>%
      dplyr::filter(
        NGS_type %in% selected_ngs,
        IMGTHLA_version %in% input$imgthla
      ) %>%
      dplyr::pull(Typing_package) %>%
      unique()
  } else {
    HLA_typing %>%
      dplyr::filter(
        NGS_type %in% input$sequence,
        IMGTHLA_version %in% input$imgthla
      ) %>%
      dplyr::pull(Typing_package) %>%
      unique()
  }
})

observe({
  if (!identical(source_mode(), "demo")) return()
  
  selected_ngs <- demo_selected_ngs()
  if (is.null(selected_ngs) || !nzchar(selected_ngs)) return()
  
  updateSelectInput(
    session,
    "sequence",
    choices = c(selected_ngs),
    selected = selected_ngs
  )
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
# ============================================================
tmp_dir <- tempdir()

observeEvent(input$step2, {
  if (identical(source_mode(), "demo")) {
    if (isTRUE(demo_uploaded())) {
      updateTabItems(session, "my_tabs", "hla")
    } else {
      shinyalert(
        title = "Error",
        text = "Please upload a demo dataset before entering the HLA typing step.",
        type = "error"
      )
    }
    return()
  }
  
  # FASTQ mode
  if (identical(source_mode(), "fastq")) {
    tmp_dir <- tempdir()
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
    return()
  }
  
  # BAM mode
  if (identical(source_mode(), "bam")) {
    tmp_dir <- tempdir()
    bam_file <- file.path(tmp_dir, "BAM", "sample.bam")
    
    if (file.exists(bam_file)) {
      updateTabItems(session, "my_tabs", "hla")
    } else {
      shinyalert(
        title = "Error",
        text = "Please upload a valid BAM file before the HLA typing step.",
        type = "error"
      )
    }
    return()
  }
})


