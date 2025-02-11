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
library(igvShiny)
library(shinyjs)
library(dplyr)
library(GenomicAlignments)
library(rtracklayer)
library(R.utils)
library(shinyalert)
library(readxl)
library(parallel)
library(openxlsx)
library(sortable)
library(summaryBox)
library(shinyWidgets)
library(knitr)
library(kableExtra)
library(formattable) 

## Upload file size limit
options(shiny.maxRequestSize = 5000 * 1024^2)

## Cascade filter
HLA_typing <- read_excel("/root/shiny/HLA_package_list.xlsx")

## System CPU
num_cores <- detectCores() - 1

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
      menuItem("Home", tabName = "home",icon = icon("house")),
      menuItem("Data Upload", tabName = "upload",icon = icon("upload")),
      menuItem("HLA Typing", tabName = "hla",icon = icon("play"))
    )
  ),
  bs4DashBody(
    shinyjs::useShinyjs(),
    autoWaiter("pivot_table",html = spin_loader(), color=NULL),
    useWaiter(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$style(
        HTML("
        .dataTables_wrapper .dataTables_paginate .paginate_button {padding: 0;margin: 0;}
        table.dataTable {width: 100% !important;}
        .dataTables_wrapper {width: 100%;overflow-x: auto;}
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
        .pvtTable {width: 100% !important;table-layout: fixed;}
        .pvtTable td, .pvtTable th {overflow: hidden;text-overflow: ellipsis;white-space: nowrap;}
      "),
        HTML("
        .radio-inline {
          margin-right: 20px;
        }
      ")
      ),
      tags$script(HTML("
      $(window).on('resize', function(){
        Shiny.onInputChange('window_resized', Math.random());
      });
      Shiny.addCustomMessageHandler('igv_resized', function(message) {
        $('#igvShiny_I').height($(window).height() * 0.6);
        $('#igvShiny_II').height($(window).height() * 0.6);
      });
    "))
    ),
    tabItems(
      tabItem(tabName = "home",source("UI/1_home.R", local = TRUE)$value),
      tabItem(tabName = "upload", source("UI/2_upload.R", local = TRUE)$value),
      tabItem(tabName = "hla", source("UI/3_hla_typing.R", local = TRUE)$value)
    )
  )
)

## Server
server <- function(input, output, session) {
  ## home page
  source("Server/1_home_html.R",local = TRUE)$vale
  
  ## Demo upload
  source("Server/2_Demo_upload.R",local = TRUE)$vale
  
  ## upload
  source("Server/2_upload_file.R",local = TRUE)$value
  
  ## QC
  source("Server/2_QC.R",local = TRUE)$vale
  
  ## cascade
  source("Server/3_cascade_filter.R",local = TRUE)$value
  
  ## package
  source("Server/3_package.R",local = TRUE)$value
  
  ##pivottabler
  source("Server/3_pivottabler.R",local = TRUE)$value
  
  ## BAM & IGV 
  source("Server/3_igvshiny.R",local = TRUE)$value

  ## valuebox
  source("Server/3_valuebox.R",local = TRUE)$value
}

shinyApp(ui, server)
