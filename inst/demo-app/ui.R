library(shiny)
library(shinyjs)
library(shinyAce)
library(dygraphs)
library(linpk)

fluidPage(
  useShinyjs(),

  #titlePanel("Simulate a Concentration-Time Profile"),

  tabsetPanel(id="toptab",
    #tabPanel("Plot", plotOutput("plot")),
    tabPanel("Plot", dygraphOutput("plot")),
    tabPanel("Table", DT::dataTableOutput("secondary")),
    tabPanel("Download",
      DT::dataTableOutput("pkprofile", width="400px"),
      downloadButton("download_btn")
      ),
    tabPanel("Code", aceEditor("code", value="", theme="xcode", mode="r",
                               readOnly=TRUE, highlightActiveLine=FALSE))
    ),


  hr(),

  wellPanel(
    fluidRow(
      column(1,
        h3("Scales"),
        numericInput("timerange", "Time range (max)", value=10, min=0, step=1),
        numericInput("ntimepoints", "Time points (N)", value=1000, min=0, step=1),
        selectInput("timeu", "Time units",
          choices=c("hours", "days", "weeks"), selected="days"),
        selectInput("doseu", "Dose units",
          choices=c("ng"="ng", "\U03BCg"="ug", "mg"="mg"), selected="mg"),
        selectInput("concu", "Concentration units",
          choices=c("ng/mL"="ng/mL", "\U03BCg/mL"="ug/mL", "mg/mL"="mg/mL"), selected="ug/mL")
        ),
      column(3,
        h3("Parameters"),
        helpText("Enter a starting value, then use the sliders to vary by up to \U00B1 100%"),
        #h4("Central Compartment"),
        fluidRow(
          column(12, checkboxInput("haska", "Absorption compartment?", TRUE))
          ),
        shinyjs::hidden(
          div(id="showka",
            fluidRow(
              column(3, numericInput("kastart", "Ka (1/d)", 10, step=1, min=0)),
              column(9, sliderInput("ka", "", min=0, max=2*10, value=10, tick=FALSE))
              )
            )
          ),
        fluidRow(
          column(3, numericInput("clstart", "CL (L/d)", 10, step=1, min=0)),
          column(9, sliderInput("cl", "", min=0, max=2*10, value=10, tick=FALSE))
          ),
        fluidRow(
          column(3, numericInput("vcstart", "VC (L)", 5, step=1, min=0)),
          column(9, sliderInput("vc", "", min=0, max=2*5, value=5, tick=FALSE))
          ),
        #h4("Peripheral Compartment(s)"),
        tags$label("Number of peripheral compartments"),
        fluidRow(
          column(3, numericInput("nperiph", NULL, value=0, min=0, max=2))
          ),
        shinyjs::hidden(
          div(id="periph1",
            fluidRow(
              column(3, numericInput("qstart", "Q (L/d)", 2, step=1, min=0)),
              column(9, sliderInput("q", "", min=0, max=2*2.5, value=2.5, tick=FALSE))
              ),
            fluidRow(
              column(3, numericInput("vpstart", "VP (L)", 10, step=1, min=0)),
              column(9, sliderInput("vp", "", min=0, max=2*10, value=10, tick=FALSE))
              )
            )
          ),
        shinyjs::hidden(
          div(id="periph2",
            fluidRow(
              column(3, numericInput("q2start", "Q2 (L/d)", 2, step=1, min=0)),
              column(9, sliderInput("q2", "", min=0, max=2*2.5, value=2.5, tick=FALSE))
              ),
            fluidRow(
              column(3, numericInput("vp2start", "VP2 (L)", 10, step=1, min=0)),
              column(9, sliderInput("vp2", "", min=0, max=2*10, value=10, tick=FALSE))
              )
            )
          )
        ),
      column(7,
        h3("Dosing Information"),
        tags$table(id="dose_table",
          tags$tr(
            tags$th(tags$label(id="ss_lab",    HTML("Steady state?"))),
            tags$th(tags$label(id="time_lab",  HTML("Time"))),
            tags$th(tags$label(id="amt_lab",   HTML("Amount"))),
            tags$th(tags$label(id="cmt_lab",   HTML("Compartment"))),
            tags$th(tags$label(id="ndose_lab", HTML("Number<br/>of doses"))),
            tags$th(tags$label(id="ii_lab",    HTML("Interdose<br/>interval"))),
            tags$th(tags$label(id="dur_lab",   HTML("Infusion<br/>duration"))),
            tags$th(tags$label(id="lag_lab",   HTML("Lag time"))),
            tags$th(tags$label(id="f_lab",     HTML("Bioavailable<br/>fraction")))
          #  ),
          #tags$tr(
          #  tags$td(checkboxInput("ss",    NULL)),
          #  tags$td(numericInput("t.dose", NULL, 0,         step=0.5,  min=0)),
          #  tags$td(numericInput("amt",    NULL, 100,       step=1,    min=0)),
          #  tags$td(selectInput("cmt",     NULL, choices=c("Default"=0), selected=0)),
          #  tags$td(numericInput("ndose",  NULL, 4,         step=1,    min=1)),
          #  tags$td(numericInput("ii",     NULL, 1,         step=1,    min=0)),
          #  tags$td(numericInput("dur",    NULL, 0,         step=0.1,  min=0)),
          #  tags$td(numericInput("lag",    NULL, 0,         step=0.1,  min=0)),
          #  tags$td(numericInput("f",      NULL, 1,         step=0.01, min=0))
            )
          ),
        actionButton("dose_add_btn", "Add"),
        actionButton("dose_del_btn", "Delete Last"),
        tags$style(type='text/css', "#dose_table { table-layout: fixed; }"),
        tags$style(type='text/css', "#dose_table tr th { padding: 0px 7px; vertical-align: bottom; min-width: 112px}"),
        tags$style(type='text/css', "#dose_table tr td { padding: 0px 7px; vertical-align: middle; min-width: 112px}"),
        tags$style(type='text/css', "#dose_table .selectize-control { margin-bottom: -5px; }"),
        tags$style(type='text/css', "#dose_table .form-group { margin: 0px; }"),
        tags$style(type='text/css', "#dose_table .checkbox { text-align: center; }"),
        tags$style(type='text/css', "#ss_lab { text-align: center; }")
        )
      )
    )
  )

# vim: ts=2 sw=2

