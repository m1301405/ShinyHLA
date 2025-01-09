output$home_html <- renderUI({
  tags$iframe(
    src = "welcome.html",
    style = "width:100%; height:100vh; border:none;",
    frameborder = "0"
  )
})