library(shiny)
library(shinydashboard)
library(bs4Dash)
library(DT)
library(rjson)
library(pivottabler)
library(readr)
library(htmlwidgets)
library(waiter)
library(rvest)
library(shinyjs)
library(dplyr)
library(R.utils)
library(shinyalert)
library(parallel)
library(openxlsx)
library(sortable)
library(summaryBox)
library(shinyWidgets)
library(knitr)
library(kableExtra)
library(formattable) 
library(reactable) 
library(htmltools)
library(ggplot2)

source("/root/shiny/R/app_globals.R")
source("/root/shiny/R/app_helpers.R")

## UI
ui <- bs4DashPage(
  help = T,
  dashboardHeader(
    title = "",
    center = tags$h2(
      style = "font-family:'Comic Sans MS', sans-serif;",
      "ShinyHLA: An R/Shiny-Based Web Application for HLA Typing of NGS Data"
    ),
    status = "light"
  ),
  bs4DashSidebar(
    skin = "light",
    status = "primary",
    elevation = 3,
    collapsed = TRUE,
    sidebarMenu(
      id = "my_tabs",
      menuItem("Home", tabName = "home", icon = icon("house")),
      menuItem("Documentation", tabName = "doc", icon = icon("file-lines")),
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("HLA Typing", tabName = "hla", icon = icon("play"))
    )
  ),
  bs4DashBody(
    shinyjs::useShinyjs(),
    autoWaiter("pivot_table", html = spin_loader(), color = NULL),
    useWaiter(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$style(
        HTML("
        .dataTables_wrapper .dataTables_paginate .paginate_button {padding: 0; margin: 0;}
        table.dataTable {width: 100% !important;}
        .dataTables_wrapper {width: 100%;}
      "),
        HTML("
        Shiny.addCustomMessageHandler('keepProgressBar', function(message) {
          var progressBar = document.getElementById('fastq_gz_upload_progress');
          if (progressBar) {
            progressBar.getElementsByClassName('progress-bar')[0].style.width = '100%';
            progressBar.getElementsByClassName('progress-bar')[0].textContent = '100%';
          }
        });
      "),
        HTML("
        .value-box .value {font-size: 30px;}
        .value-box .text {font-size: 18px;}
      "),
        HTML("
        .pvtTable {width: 100% !important; table-layout: fixed;}
        .pvtTable td, .pvtTable th {overflow: hidden; text-overflow: ellipsis; white-space: nowrap;}
      "),
        HTML("
        .radio-inline {
          margin-right: 20px;
        }
      "),
        HTML("
        .summary-box {
          min-height: 120px !important;
          display: flex;
          flex-direction: column;
          justify-content: center;
        }
        #hla_typing_table {
          width: 100% !important;
          margin: 0 auto;
        }
        .dataTables_filter input {
          border-radius: 20px;
          padding: 5px 15px;
          border: 1px solid #ddd;
        }
      "),
        HTML("
        /* Result tab horizontal scroll */
        .hla-result-scroll{
          width: 100%;
          overflow-x: auto;
          overflow-y: visible;
          padding-bottom: 6px;
        }
        .hla-result-inner{
          min-width: 1100px;
        }
        .hla-result-inner .row{
          margin-left: 0 !important;
          margin-right: 0 !important;
        }
        #hla_box .card-body{
          overflow: visible !important;
        }
        #hla_box .dataTables_wrapper{
          overflow: visible !important;
        }
        #summary_boxes .summary-scroll{
          width: 100%;
          overflow-x: auto;
          overflow-y: hidden;
        }
        
        #summary_boxes .summary-strip{
          display: flex;
          flex-wrap: nowrap;
          width: 100%;
        }
        
        #summary_boxes .summary-item{
          flex: 1 1 25%;
          min-width: 260px;
          padding: 0 8px;
          box-sizing: border-box;
        }
        
        #summary_boxes .summary-box{
          width: 100% !important;
          min-height: 110px !important;
          margin-bottom: 0 !important;
        }
      "),
        HTML("
        /* ================= PharmGKB embedded styles ================= */
        .summary-card {
          background: white;
          border: 1px solid #e5e5e5;
          border-radius: 12px;
          padding: 18px 20px;
          margin-bottom: 18px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.04);
        }
        
        .summary-badge {
          display: inline-block;
          padding: 8px 16px;
          background: #fdecea;
          color: #b71c1c;
          border-radius: 18px;
          font-weight: 700;
          font-size: 14px;
        }
        
        .mini-card {
          background: #fafafa;
          border: 1px solid #eee;
          border-radius: 10px;
          padding: 14px;
          min-width: 180px;
        }
        
        .mini-title {
          font-size: 12px;
          color: #777;
          text-transform: uppercase;
        }
        
        .mini-value {
          font-size: 26px;
          font-weight: 700;
        }
        
        .overview-flex {
          display: flex;
          flex-wrap: nowrap;
          gap: 28px;
          align-items: flex-start;
        }
        
        .overview-left {
          flex: 0 0 58%;
          min-width: 0;
        }
        
        .overview-right {
          flex: 0 0 38%;
          min-width: 340px;
        }
        
        .overview-metrics {
          display: grid;
          grid-template-columns: repeat(3, minmax(0, 1fr));
          gap: 14px;
          margin-bottom: 18px;
        }
        
        @media (max-width: 1100px) {
          .overview-flex {
            flex-wrap: wrap;
          }
          
          .overview-left,
          .overview-right {
            flex: 0 0 100%;
            min-width: 100%;
          }
          
          .overview-right {
            margin-top: 8px;
          }
        }
        
        @media (max-width: 768px) {
          .overview-metrics {
            grid-template-columns: 1fr;
          }
        }
      "),
        HTML("
        #pharmgkb_box {
          width: 100% !important;
        }
      
        #pharmgkb_box .card-body {
          width: 100% !important;
        }
      
        #pharmgkb_box .tab-content {
          width: 100% !important;
        }
      
        #pharmgkb_box .tab-pane {
          width: 100% !important;
          min-width: 100% !important;
        }
      
        #pharmgkb_box .summary-card {
          width: 100% !important;
          box-sizing: border-box;
        }
      ")
      ),
      tags$script(HTML("
      $(window).on('resize', function(){
        Shiny.onInputChange('window_resized', Math.random());
      });
    
      Shiny.addCustomMessageHandler('resetPharmgkbButtons', function(message) {
        document.querySelectorAll('.pharmgkb-search-btn').forEach(function(btn) {
          btn.disabled = false;
          btn.innerHTML = 'Search';
        });
      });
    "))
    ),
    tabItems(
      tabItem(tabName = "home", source("UI/1_home.R", local = TRUE)$value),
      tabItem(tabName = "doc", source("UI/2_doc.R", local = TRUE)$value),
      tabItem(tabName = "upload", source("UI/3_upload.R", local = TRUE)$value),
      tabItem(tabName = "hla", source("UI/4_hla_typing.R", local = TRUE)$value)
    )
  )
)

## Server
server <- function(input, output, session) {
  # --- global reactive states (single source of truth) ---
  current_job_id      <- reactiveVal(NULL)
  current_tmp_zip     <- reactiveVal(NULL)
  current_history_zip <- reactiveVal(NULL)
  current_mode        <- reactiveVal("run")
  package_value   <- reactiveVal("")
  current_zippath <- reactiveVal(NULL)
  current_prefix  <- reactiveVal(NULL)
  HISTORY_ROOT <- "/root/shiny/History"
  dir.create(HISTORY_ROOT, recursive = TRUE, showWarnings = FALSE)
  
  ## Demo upload
  source("Server/3_Demo_cascade.R",local = TRUE)$value
  source("Server/3_Demo_upload.R",local = TRUE)$value
  source("Server/3_Demo_pivotable.R",local = TRUE)$value
  
  ## upload
  source("Server/3_upload_file.R",local = TRUE)$value
  
  ## QC
  source("Server/3_QC.R",local = TRUE)$vale
  
  ## cascade
  source("Server/4_cascade_filter.R",local = TRUE)$value
  
  ## package
  source("Server/4_package.R",local = TRUE)$value
  
  ## pivottabler
  source("Server/4_pivottabler.R",local = TRUE)$value
  
  ## pharmgkb
  source("Server/4_pharmgkb.R", local = TRUE)$value
  
  ## valuebox
  source("Server/4_valuebox.R",local = TRUE)$value
}

shinyApp(ui, server)
