library(shiny)
library(shinyjs)
library(dygraphs)
library(linpk)

function(input, output, session) {

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
    shinyjs::hide("iicontainer")
    if (input$ss || input$ndose > 1) {
      shinyjs::show("iicontainer")
    }
  })

  observe({
    shinyjs::hide("showka")
    if (input$haska) {
      shinyjs::show("showka")
    }
  })

  simulated <- reactive({
    validate(
      need(as.numeric(input$cl) > 0, "Enter a positive number for CL"),
      need(as.numeric(input$vc) > 0, "Enter a positive number for VC"),
      need(input$nperiph < 1 | as.numeric(input$q) > 0, "Enter a positive number for Q"),
      need(input$nperiph < 1 | as.numeric(input$vp) > 0, "Enter a positive number for VP"),
      need(input$nperiph < 2 | as.numeric(input$q2) > 0, "Enter a positive number for Q2"),
      need(input$nperiph < 2 | as.numeric(input$vp2) > 0, "Enter a positive number for VP2"),
      need(input$haska == FALSE | as.numeric(input$ka) > 0, "Enter a positive number for Ka"),
      need(as.numeric(input$lag) >= 0, "Enter a non-negative number for lag time"),
      need(as.numeric(input$f) > 0, "Enter a positive number for lag bioavailable fraction")
      )

    t.obs <- seq(0, input$timerange, length.out=1000)
    cl <- as.numeric(input$cl)
    vc <- as.numeric(input$vc)
    q <- as.numeric(c(input$q, input$q2))[seq_len(input$nperiph)]
    vp <- as.numeric(c(input$vp, input$vp2))[seq_len(input$nperiph)]
    t.dose <- as.numeric(input$t.dose)
    addl <- as.numeric(input$ndose) - 1
    ss <- as.numeric(input$ss)
    ii <- as.numeric(input$ii)
    amt <- as.numeric(input$amt)/switch(input$doseu, "ng"=10^6, "mg"=1, 10^3)
    lag <- as.numeric(input$lag)
    f <- as.numeric(input$f)
    ka <- 0
    if (input$haska) {
      ka <- as.numeric(input$ka)
    }
    dur <- as.numeric(input$dur)

    y <- pkprofile(t.obs, cl=cl, vc=vc, q=q, vp=vp, ka=ka, dose=list(t.dose=t.dose, amt=amt, addl=addl, ii=ii, ss=ss, ka=ka, dur=dur, lag=lag, f=f))
    data.frame(time=t.obs, conc=as.numeric(y)*switch(input$concu, "ng/mL"=10^3, "mg/mL"=10^-3, 1))
  })

  output$dygraph <- renderDygraph({
    dygraph(simulated(), main="PK Concentration-Time Profile") %>%
      dySeries("conc", label="Model 1") %>%
      dyOptions(drawGrid=TRUE, includeZero=TRUE) %>%
      dyAxis("y", sprintf("Concentration (%s)", input$concu)) %>%
      dyAxis("x", sprintf("Time (%s)", input$timeu))
  })

}

# vim: ts=2 sw=2

