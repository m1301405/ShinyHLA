package_value <- reactiveVal()

observeEvent(input$hla_typing_button, {
  waiter <- Waiter$new(
    id = "summary_boxes",
    html = spin_loader(),  
    color = transparent(.5)
  )
  waiter$show()
  
  NGS_value <- isolate(input$sequence)
  imgthla_value <- isolate(input$imgthla)
  
  package_value(isolate(input$package))
  
  render_summary_boxes_ui <- function(ngs, version, tool, jid) {
    div(
      class = "summary-scroll",
      div(
        class = "summary-strip",
        div(class = "summary-item",
            summaryBox("NGS TYPE", ngs, width = 12, icon = icon("dna"), style = "info")
        ),
        div(class = "summary-item",
            summaryBox("VERSION", version, width = 12, icon = icon("code-branch"), style = "success")
        ),
        div(class = "summary-item",
            summaryBox("TYPING TOOL", tool, width = 12, icon = icon("wrench"), style = "danger")
        ),
        div(class = "summary-item",
            summaryBox("JOB ID", jid, width = 12, icon = icon("id-card"), style = "primary")
        )
      )
    )
  }
  
  observeEvent(input$hla_typing_button, {
    package_value(isolate(input$package))  # 讓 UI/conditionalPanel 依賴一致
    
    output$summary_boxes <- renderUI({
      render_summary_boxes_ui(
        ngs     = isolate(input$sequence),
        version = isolate(input$imgthla),
        tool    = isolate(input$package),
        jid     = current_job_id()
      )
    })
  })
  
  observeEvent(package_value(), {
    updateTextInput(session, "package_value_ui", value = package_value())
  })
  
  waiter$hide()
})

observeEvent(package_value(), {
  updateTextInput(session, "package_value_ui", value = package_value())
})