output$home_html <- renderUI({
  mtime <- as.integer(file.info("www/welcome.html")$mtime)
  
  tags$iframe(
    src = sprintf("welcome.html?v=%s", mtime),
    style = "width:100%; height:100vh; border:none;",
    frameborder = "0"
  )
})