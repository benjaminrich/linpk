library(shiny)
library(shinyjs)
library(shinyAce)
library(dygraphs)
library(linpk)

function(input, output, session) {

  v <- reactiveValues(num_dose=0, next.dose.t=0)

  observe({
    updateSliderInput(session, "ka",
      value=input$kastart, min=0, max=2*input$kastart, step=input$kastart/50)
  })

  observe({
    updateSliderInput(session, "cl",
      value=input$clstart, min=0, max=2*input$clstart, step=input$clstart/50)
  })

  observe({
    updateSliderInput(session, "vc",
      value=input$vcstart, min=0, max=2*input$vcstart, step=input$vcstart/50)
  })

  observe({
    updateSliderInput(session, "q",
      value=input$qstart, min=0, max=2*input$qstart, step=input$qstart/50)
  })

  observe({
    updateSliderInput(session, "vp",
      value=input$vpstart, min=0, max=2*input$vpstart, step=input$vpstart/50)
  })

  observe({
    updateSliderInput(session, "q2",
      value=input$q2start, min=0, max=2*input$q2start, step=input$q2start/50)
  })

  observe({
    updateSliderInput(session, "vp2",
      value=input$vp2start, min=0, max=2*input$vp2start, step=input$vp2start/50)
  })

  observe({
    shinyjs::hide("periph1")
    shinyjs::hide("periph2")
    if (input$nperiph >=1) {
      shinyjs::show("periph1")
    }
    if (input$nperiph >=2) {
      shinyjs::show("periph2")
    }
  })

  observe({
    shinyjs::hide("showka")
    if (input$haska) {
      shinyjs::show("showka")
    }
  })

  observe({
    short <- switch(input$timeu, "hours"="h", "days"="d", "weeks"="w")
    updateNumericInput(session, "kastart", label=sprintf("Ka (1/%s)", short))
    updateNumericInput(session, "clstart", label=sprintf("CL (L/%s)", short))
    updateNumericInput(session, "qstart", label=sprintf("Q (L/%s)", short))
    updateNumericInput(session, "q2start", label=sprintf("Q2 (L/%s)", short))
  })

  compartments <- reactive({
    res <- c("Default"=0, "Central"=1)
    if (input$nperiph == 1) {
      res <- c(res, "Peripheral"=2)
    } else if (input$nperiph == 2) {
      res <- c(res, "Periph. 1"=2, "Periph. 2"=3)
    }
    if (input$haska) {
      res <- c(res, "Absorption"=max(res)+1)
    }
    res
  })

  observeEvent(input$dose_add_btn, ignoreNULL=FALSE, {
    i <- v$num_dose + 1
    t <- v$next.dose.t 
    insertUI(
      selector="#dose_table", where="beforeEnd", immediate=TRUE,
      tags$tr(id=paste0("dose_row_", i),
        tags$td(checkboxInput(paste0("ss_", i),    NULL)),
        tags$td(numericInput(paste0("t.dose_", i), NULL, t,   step=0.5,  min=0)),
        tags$td(numericInput(paste0("amt_", i),    NULL, 100, step=1,    min=0)),
        tags$td(selectInput(paste0("cmt_", i),     NULL, selected=0, choices=compartments())),
        tags$td(numericInput(paste0("ndose_", i),  NULL, 1,   step=1,    min=1)),
        tags$td(numericInput(paste0("ii_", i),     NULL, 1,   step=1,    min=0)),
        tags$td(numericInput(paste0("dur_", i),    NULL, 0,   step=0.1,  min=0)),
        tags$td(numericInput(paste0("lag_", i),    NULL, 0,   step=0.1,  min=0)),
        tags$td(numericInput(paste0("f_", i),      NULL, 1,   step=0.01, min=0))
        )
      )
    v$num_dose <- i
  })

  observeEvent(input$dose_del_btn, {
    if (v$num_dose > 1) {
      i <- v$num_dose
      removeUI(selector=paste0("#dose_row_", i))
      v$num_dose <- v$num_dose - 1
    }
  })

  observe({
    toggleState("dose_del_btn", condition=(v$num_dose > 1))
  })

  observe({
    for (i in seq_len(v$num_dose)) {
      updateSelectInput(session, paste0("cmt_", i), choices=compartments())
    }
  })


  sim <- reactive({
    t.obs <- seq(0, input$timerange, length.out=input$ntimepoints)
    cl    <- as.numeric(input$cl)
    vc    <- as.numeric(input$vc)
    q     <- as.numeric(c(input$q, input$q2))[seq_len(input$nperiph)]
    vp    <- as.numeric(c(input$vp, input$vp2))[seq_len(input$nperiph)]
    ka    <- 0
    if (input$haska) {
      ka <- as.numeric(input$ka)
    }

    validate(
      need(cl > 0, "Enter a positive number for CL"),
      need(vc > 0, "Enter a positive number for VC"),
      need(input$nperiph < 1 | q[1] > 0, "Enter a positive number for Q"),
      need(input$nperiph < 1 | vp[1] > 0, "Enter a positive number for VP"),
      need(input$nperiph < 2 | q[2] > 0, "Enter a positive number for Q2"),
      need(input$nperiph < 2 | vp[2] > 0, "Enter a positive number for VP2"),
      need(input$haska == FALSE | ka > 0, "Enter a positive number for Ka")
      )

    validate(need(v$num_dose > 0, "No dose information"))

    dose <- data.frame(
      ss     = rep(FALSE, v$num_dose),
      t.dose = 0,
      amt    = 0,
      cmt    = 0,
      addl   = 0,
      ii     = 0,
      dur    = 0,
      lag    = 0,
      f      = 0)

    try({
      for (i in seq_len(v$num_dose)) {
        dose$ss[i]     <- as.logical(input[[paste0("ss_", i)]])
        dose$t.dose[i] <- as.numeric(input[[paste0("t.dose_", i)]])
        dose$amt[i]    <- as.numeric(input[[paste0("amt_", i)]])
        dose$cmt[i]    <- as.numeric(input[[paste0("cmt_", i)]])
        dose$addl[i]   <- as.numeric(input[[paste0("ndose_", i)]]) - 1
        dose$ii[i]     <- as.numeric(input[[paste0("ii_", i)]])
        dose$dur[i]    <- as.numeric(input[[paste0("dur_", i)]])
        dose$lag[i]    <- as.numeric(input[[paste0("lag_", i)]])
        dose$f[i]      <- as.numeric(input[[paste0("f_", i)]])
      }
    }, TRUE)

    v$next.dose.t <- dose$t.dose[i] + max(dose$ii[i], 1)

    dose$amt <- dose$amt/switch(input$doseu, "ng"=10^6, "mg"=1, 10^3)
    dose$amt <- dose$amt*switch(input$concu, "ng/mL"=10^3, "mg/mL"=10^-3, 1)

    validate(
      need(all(dose$lag >= 0), "Enter a non-negative number for lag time"),
      need(all(dose$f > 0), "Enter a positive number for lag bioavailable fraction")
      )

    tryCatch({
      pkprofile(t.obs, cl=cl, vc=vc, q=q, vp=vp, ka=ka, dose=dose)
    },
    error = function(err) {
      err$message
    })
  })

  sim.df <- reactive({
    y <- sim()
    if (is.character(y)) {
      return(y)
    }
    as.data.frame(y)
  })

  #output$plot <- renderPlot({
  #  validate(need(!is.null(sim.df()), "Invalid data"))
  #  plot(conc ~ time, data=sim.df(), type="l", col="blue", lwd=1.1,
  #    ylim=c(0, max(sim())),
  #    main="PK Concentration-Time Profile",
  #    ylab=sprintf("Concentration (%s)", gsub("ug", "\U03BCg", input$concu)),
  #    xlab=sprintf("Time (%s)", input$timeu))
  #})

  output$plot <- renderDygraph({
    df <- sim.df()
    if (is.character(df)) {
      stop(df)
    }
    validate(need(!is.null(df), "Invalid data"))
    dygraph(df, main="PK Concentration-Time Profile") %>%
      dySeries("conc", label="Model 1") %>%
      dyOptions(drawGrid=TRUE, includeZero=TRUE) %>%
      dyAxis("y", sprintf("Concentration (%s)", gsub("ug", "&mu;g", input$concu))) %>%
      dyAxis("x", sprintf("Time (%s)", input$timeu))
  })

  simtab <- reactive({
    y <- sim()
    if (is.character(y)) {
      return(y)
    }
    as.data.frame(secondary(y))
  })

  output$secondary <- DT::renderDataTable({
    tab <- simtab()
    if (is.character(tab)) {
      stop(tab)
    }
    format(tab)
  }, options=list(dom="t"))

  output$pkprofile <- DT::renderDataTable({
    df <- sim.df()
    if (is.character(df)) {
      stop(df)
    }
    format(df)
  }, options=list(dom="tip"))

  output$download_btn <- downloadHandler(
    filename="pkprofile.csv",
    content=function(file) {
      tryCatch({
        write.csv(format(sim.df()), file, row.names=FALSE)
      },
      error = function(err) {
        stop(err$message)
      })
    })

  observe({
    input$toptab  # Create dependency
    args <- character(0)
    args <- c(args, sprintf("t.obs = seq(0, %s, length.out=%s)", input$timerange, input$ntimepoints))
    args <- c(args, sprintf("cl = %s", input$cl))
    args <- c(args, sprintf("vc = %s", input$vc))
    if (input$nperiph == 1) {
      args <- c(args, sprintf("q = %s", input$q))
      args <- c(args, sprintf("vp = %s", input$vp))
    } else if (input$nperiph > 1) {
      args <- c(args, sprintf("q = c(%s, %s)", input$q, input$q2))
      args <- c(args, sprintf("vp = c(%s, %s)", input$vp, input$vp2))
    }
    if (input$haska) {
      args <- c(args, sprintf("ka = %s", input$ka))
    }
    args <- c(args, "dose = dose")

    f <- function(x) {
      if (length(x) == 0) {
        "NULL"
      } else if (length(x) == 1) {
        as.character(x)
      } else {
        sprintf("c(%s)", paste(x, collapse=", "))
      }
    }

    doseargs <- character(0)
    y <- sim()
    dose <- attr(y, "dose")
    doseargs <- c(doseargs, sprintf("t.dose = %s", f(dose$t.dose)))
    doseargs <- c(doseargs, sprintf("amt = %s", f(dose$amt)))
    if (!is.null(dose$cmt) && any(dose$cmt != 0)) {
      doseargs <- c(doseargs, sprintf("cmt = %s", f(dose$cmt)))
    }
    if (!is.null(dose$addl) && any(dose$addl != 0)) {
      doseargs <- c(doseargs, sprintf("addl = %s", f(dose$addl)))
    }
    if (!is.null(dose$ii) && any(dose$addl != 0 | dose$ss)) {
      doseargs <- c(doseargs, sprintf("ii = %s", f(dose$ii)))
    }
    if (!is.null(dose$dur) && any(dose$dur != 0)) {
      doseargs <- c(doseargs, sprintf("dur = %s", f(dose$dur)))
    }
    if (!is.null(dose$lag) && any(dose$lag != 0)) {
      doseargs <- c(doseargs, sprintf("lag = %s", f(dose$lag)))
    }
    if (!is.null(dose$f) && any(dose$f != 1)) {
      doseargs <- c(doseargs, sprintf("f = %s", f(dose$f)))
    }
    if (!is.null(dose$ss) && any(dose$ss)) {
      doseargs <- c(doseargs, sprintf("ss = %s", f(dose$ss)))
    }

    codestr <- "require(linpk)"
    codestr <- c(codestr, sprintf("dose <- data.frame(%s)", paste(doseargs, collapse=", ")))
    codestr <- c(codestr, sprintf("y <- pkprofile(%s)", paste(args, collapse=", ")))
    codestr <- paste(codestr, collapse="\n")
    updateAceEditor(session, "code", codestr)
  })
}

# vim: ts=2 sw=2

