# 定義全局變量 package_value，使用 reactiveVal 來儲存
# 注意：package_value 與 current_job_id 應在 app.R 或 server 頂層定義
package_value <- reactiveVal()

observeEvent(input$hla_typing_button, {
  # 創建等待動畫
  waiter <- Waiter$new(
    id = "summary_boxes",
    html = spin_loader(),  # 使用 spin_loader() 動畫
    color = transparent(.5)
  )
  # 顯示等待動畫
  waiter$show()
  
  # 取得輸入值
  NGS_value <- isolate(input$sequence)
  imgthla_value <- isolate(input$imgthla)
  
  # 更新全局變量 package_value
  package_value(isolate(input$package))
  
  # 設定 summaryBox UI
  # 修改重點：將 width 從 4 改為 3，並增加第四個 JobID 欄位
  output$summary_boxes <- renderUI({
    fluidRow(
      # 第一個箱子：NGS Type
      column(width = 3, summaryBox("NGS Type", NGS_value, width = 12, icon = icon("dna"), style = "info")),
      
      # 第二個箱子：IPD-IMGT/HLA Version
      column(width = 3, summaryBox("IPD-IMGT/HLA Version", imgthla_value, width = 12, icon = icon("code-branch"), style = "success")),
      
      # 第三個箱子：Typing Tool
      column(width = 3, summaryBox("Typing Tool", package_value(), width = 12, icon = icon("box"), style = "danger")),
      
      # 第四個箱子：JobID (新增)
      # 寬度設為 3，內部 summaryBox 寬度設為 12 撐滿
      column(width = 3, summaryBox("JobID", current_job_id(), width = 12, icon = icon("fingerprint"), style = "primary"))
    )
  })
  
  # 隱藏等待動畫
  waiter$hide()
})

observeEvent(package_value(), {
  updateTextInput(session, "package_value_ui", value = package_value())
})