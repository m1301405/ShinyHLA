# Define a global variable `package_value` using reactiveVal for storage
package_value <- reactiveVal()

observeEvent(input$hla_typing_button, {
  # Create a loading animation
  waiter <- Waiter$new(
    id = "summary_boxes",
    html = spin_loader(),  # Use the spin_loader() animation
    color = transparent(.5)
  )
  # Display the loading animation
  waiter$show()
  
  # Retrieve input values
  NGS_value <- isolate(input$sequence)
  imgthla_value <- isolate(input$imgthla)
  
  # Update the global variable `package_value`
  package_value(isolate(input$package))
  
  # Set the UI for the summary boxes
  output$summary_boxes <- renderUI({
    fluidRow(
      column(width = 4, summaryBox("NGS Type", NGS_value, width = 12, icon = icon("dna"), style = "info")),
      column(width = 4, summaryBox("IPD-IMGT/HLA Version", imgthla_value, width = 12, icon = icon("code-branch"), style = "success")),
      column(width = 4, summaryBox("Typing Tool", package_value(), width = 12, icon = icon("box"), style = "danger"))
    )
  })
  
  # Hide the loading animation
  waiter$hide()
})

observeEvent(package_value(), {
  updateTextInput(session, "package_value_ui", value = package_value())
})
