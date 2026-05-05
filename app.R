library(shiny)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(mrgsolve)
library(minqa)
library(MASS)
library(phonTools)
library(gridExtra)

source('helpers_IFX.R')

# infList <- list(
#   data.frame(dose = 20, dt = as.POSIXct("2022-01-01 00:00:00"), lab = 1, valid = TRUE),
#   data.frame(dose = 10, dt = as.POSIXct("2022-01-01 00:00:00"), lab = 1, valid = TRUE)
# )
# conList <- data.frame(dose = 20, dt = as.POSIXct('2022-01-01 4:30:00'), lab = 1, valid = TRUE)
# ex_res <- setupModel(infList, conList, 'infliximab', wt=70, alb=4, ada=0)
#
# df <- read.csv(stdin())
# rate,amt,conc,dt,lime,duration
# 20,20,NA,"2022-01-01 00:00:00",0,1
# NA,NA,20,"2022-01-01 04:30:00",4.5,NA
# NA,1,NA,"2022-01-01 10:30:00",10.5,NA
# 
# do.call(setupModel, c(arrange_data(df), drug='vancomycinf', wt=2.9, age=10, pma=10, creat=0.4))

qweek <- sprintf('Q%sWK', 1:8)
h1 <- tags$div(textInput('targetTrough', 'Target Trough (mcg/mL)', width = '150px'), style = "display:inline-block")
h2 <- tags$div(actionButton('findTarget', 'Find'), style = "display:inline-block")

ui <- fluidPage(
  useShinyjs(),
  
  tags$div(
    style = "color:#8B0000; font-weight:bold; border:1px solid #8B0000; padding:10px; margin-bottom:15px; border-radius:5px;",
    "Disclaimer: These predictions are estimates and may not accurately reflect your patient population or individual patient circumstances. ",
    "This tool is not a substitute for clinical judgment, institutional protocols, or standards of care. ",
    "All dosing decisions remain the responsibility of the treating clinician. Use at your own risk."
  ),
  
  titlePanel("Therapeutic Drug Monitoring"),
  tabsetPanel(
    tabPanel("Application",
      fluidRow(
        column(3,
          hr(),
          fileInput('file1', 'Upload Dosing Profile (CSV)'),
          hr(),
          selectInput("drug", "Drug:",
                    c("Infliximab" = "infliximab")),
          hr(),
          radioButtons("unit", "Infusion Unit", choices=c("mg/kg", "mg"), selected = 'mg/kg', inline=TRUE),
          tabsetPanel(type="pills",
          tabPanel("Demographics and Labs",
            textInput("wt", "Weight (kg)", value = "70"),
            textInput("alb", "Albumin (g/dL)", value = "4"),
            textInput("ada", "Anti-drug antibody", value = "0")
          ),
          tabPanel("Dosing 1",
            htmlOutput('inf1')
          ),
          tabPanel("Dosing 2",
            htmlOutput('inf2')
          ),
          tabPanel("Dosing 3",
            htmlOutput('inf3')
          ),
          tabPanel("Drug Levels",
            htmlOutput('conc')
          )
          )
        ),
        column(6,
          tags$br(),
          actionButton('runmodel', 'Create Output'),
          plotOutput("responseplot"),
          align = 'center'
        ),
        column(3,
          h2("Induction Phase"),
          tableOutput('peaktroughtable'),
          hr(),
          h2("Maintenance Phase"),
          textInput("ud", "Custom Dose (mg/kg)", value = "10", width = "150px"),
          radioButtons("ufrq", "Frequency", choices = qweek, inline = FALSE),
          tableOutput('utable'),
          div(h1, h2),
          hr(),
          h4("Export Results"),
          downloadButton('downloadPlot', 'Download TDM Plot (PDF)'),
          downloadButton('downloadTable', 'Download Dosing Profile (CSV)'),
          style='margin-bottom:30px;leftborder:1px solid; padding: 10px;font-size:0.8em'
        )
      )
    ),
    tabPanel("Documentation",
             h1("Functionality"),
             p("This application takes as input patient characteristics/dosing schedule and presents both the population expected drug concentrations
               as well as a simulated set of drug concentrations that represent randomly generated individuals according to a model-estimated 
               random effects distribution.  These simulated responses give some intuition to how individuals can be expected to 
               deviate.  Additionally, if drug concentrations have been observed, we can use that information to calculate the empirical Bayes 
               estimate (EBE) for the random effects and plot the individual expected drug concentrations in order to aid therapeutic drug monitoring (TDM).  Further, we can specify two alternative dosing schedules and see population, simulated, and individual expected drug concentrations to those schedules."),
             h1("Resources"),
             HTML("<p>The models implemented in this application are all available in the scientific literature.
                  The infliximab model is published in
                  <a href='https://pubmed.ncbi.nlm.nih.gov/27739008/'>Dubinsky 2017</a> [1].  
                  The code to calculate the EBEs is adapted from the open source TDM software 
                  <a href='https://mrgsolve.org/'>mrgsolve</a>.
                  </p>"),
             h1("Bibliography"),
             tags$ol(
               tags$li("Dubinsky MC, Phan BL, Singh N, Rabizadeh S, Mould DR. Pharmacokinetic Dashboard-Recommended Dosing Is Different than Standard of Care Dosing in Infliximab-Treated Pediatric IBD Patients. AAPS J. 2017 Jan;19(1):215-222. doi: 10.1208/s12248-016-9994-y."), 
               tags$li("Xu Z, Mould DR, Hu C, Ford J, Keen M, Davis HM, et al. A population-based pharmacokinetic pooled analysis of infliximab in pediatrics. ACCP National Meeting 2012 San Diego CA")
             )
    ),
    tabPanel("View Profile",
      DT::dataTableOutput('profile1'),
      downloadButton('downloadTable1', 'Download Data (CSV)')
    )
  )
)

data2Inputs <- function(dd, myid) {
  doseLabel <- c('Bolus','Infusion','Conc')[match(substr(myid, 1, 1), c('b','i','c'))]
  nn <- nrow(dd)
  myclick <- sprintf('Shiny.onInputChange( \"delete_%s_button\" , this.id, {priority: \"event\"})', myid)
  divstyle <- "display:inline-block;vertical-align:top;height:40px"
  inputs <- vector('list', nn + 2)
#   ids <- paste0('delete_', dd[,'lab'])

  for(i in seq_len(nn)) {
    dose <- textInput(paste0('dose_', myid, '_', i), label = NULL, value = dd[i,'dose'], width = '60px')
    suppressMessages(
      dt <- shinyWidgets::airDatepickerInput(paste0('dt_', myid, '_', i), label = NULL, value = dd[i,'dt'], timepickerOpts = list(timeFormat = 'HH:mm'),
                                             timepicker = TRUE, width = '180px', update_on = 'close'
      )
    )
#     ab <- actionButton(ids[i], label = NULL, icon = icon('trash'), onclick = myclick)
#     inputs[[i + 1]] <- tags$div(tags$div(dose, style = divstyle), tags$div(dt, style = divstyle), tags$div(ab, style = divstyle))
    inputs[[i + 1]] <- tags$div(tags$div(dose, style = divstyle), tags$div(dt, style = divstyle))
  }
  adder <- sprintf('Shiny.onInputChange( \"add_%s_button\" , this.id, {priority: \"event\"})', myid)
  add_label <- 'add dose'
  if(doseLabel == 'Conc') add_label <- 'add conc'
  addit <- actionButton('addit', label = add_label, icon = icon('plus'), onclick = adder, style = 'padding:4px;font-size:80%')
  h1 <- tags$div(doseLabel, style = "display:inline-block;vertical-align:top;width:60px")
  h2 <- tags$div('Time', style = "display:inline-block;vertical-align:top;width:180px")
  #   h3 <- tags$div(actionButton('addit', label = 'add', icon = icon('plus'), onclick = adder, style = 'padding:2px;font-size:80%'), style = "display:inline-block;vertical-align:top")
  inputs[[1]] <- tags$div(h1, h2)
  inputs[[nn + 2]] <- tags$div(addit)
  inputs
}

server <- function(input, output) {
  v <- reactiveValues(dat = NULL, plot1 = NULL, plot2 = NULL, plot3 = NULL, params = NULL, sched = NULL)

  formulaText <- reactive({
    paste(input$drug)
  })

  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })

  observeEvent(input$file1, {
    if(!is.null(input$file1)) {
      inFile <- input$file1
      if(!is.null(inFile) && inFile$type %in% c('text/csv', 'text/comma-separated-values', 'text/plain')) {
        v$dat <- read.csv(inFile$datapath)
      } else {
        v$dat <- NULL
      }
    }
  })
  PKprof <- reactiveValues(p = NULL, Omega = NULL, wt = NULL, unit = 'mg/kg')
  observe({
    toggleElement(id = "wt", condition = is.null(v$dat$Weight))
    toggleElement(id = "downloadTable", condition = is.null(input$file1))
    toggleElement(id = "downloadPlot", condition = !is.null(v$params))
    toggleElement(id = "vp", condition = !is.null(v$params))
  })

  ####################    ####################    ####################    ####################    ####################    ####################
  def_time1 <- as.POSIXct("2022-01-01")
  def_time1 <- as.POSIXct(format(Sys.time(), "%Y-%m-%d"))
  # def_inf_seq <- seq(def_time1, length.out=2, by= 3600*24*7*2) # HERE test with magnifying
  # def_inf_seq1 <- def_inf_seq[c(1, 2)] # HERE test with magnifying
  def_inf_seq <- seq(def_time1, length.out=20, by= 3600*24*7*2 ) # HERE 3600*24*7*2: every 2 weeks ; 3600*12: 12 hour increments
  def_inf_seq1 <- def_inf_seq[c(1, 2, 4, 8)] # HERE 3600*24*7*2: every 2 weeks ; 3600*12: 12 hour increments

  def_inf_seq2 <- def_inf_seq[c(1, 3, 5, 7)]
  def_inf_seq3 <- def_inf_seq[c(1, 5, 9, 13)]

  values <- reactiveValues(
    infusion1dat = data.frame(dose = rep(5,4), dt = def_inf_seq1, lab = seq(4), valid = TRUE), # HERE change default
    # infusion2dat = data.frame(dose = rep(10,4), dt = def_inf_seq2, lab = seq(4), valid = TRUE),
    # infusion3dat = data.frame(dose = rep(10,4), dt = def_inf_seq3, lab = seq(4), valid = TRUE),
    infusion2dat = data.frame(dose = rep(7.5,4), dt = def_inf_seq1, lab = seq(4), valid = TRUE),
    infusion3dat = data.frame(dose = rep(10,4), dt = def_inf_seq1, lab = seq(4), valid = TRUE),

    concdat = data.frame(dose = 15, dt = def_time1 + 3600*(24*7*2 - 1), lab = 1, valid = TRUE), # HERE; 1 hr before next dose for obs. conc time
    init = FALSE
  )

  factory_add <- function(dat, myid) {
    button <- sprintf('add_%s_button', myid)
    observeEvent(input[[button]], {
      if(!values$init) return()
      ix <- which(values[[dat]][,'valid'])
      if(is.na(ix[1])) {
        next_time <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:00"))
      } else {
        next_time <- as.POSIXct(format(max(values[[dat]][ix,'dt']), "%Y-%m-%d %H:00")) + 3600 * 12
      }
      values[[dat]] <- rbind(values[[dat]], data.frame(dose = 0, dt = next_time, lab = nrow(values[[dat]]) + 1, valid = TRUE))
    })
  }
  factory_del <- function(dat, myid) {
#     button <- sprintf('delete_%s_button', myid)
#     observeEvent(input[[button]], {
#       if(!values$init) return()
#       selectedId <- as.numeric(strsplit(input[[button]], "_")[[1]][2])
#       selectedRow <- match(selectedId, values[[dat]][,'lab'])
#       values[[dat]][selectedRow, 'valid'] <- FALSE
#     })
  }
  factory_dose <- function(dat, myid) {
    doseid <- sprintf('dose_%s_', myid)
    observeEvent({jj <- nrow(values[[dat]]); lapply(paste0(doseid, seq(jj)), function(i) input[[i]])}, {
      if(!values$init) return()
      labix <- values[[dat]][,'lab']
      curr_dose <- values[[dat]][,'dose']
      dose_vals <- character(length(labix))
      for(i in seq_along(dose_vals)) {
        tmp <- input[[paste0(doseid, i)]]
        if(is.null(tmp)) tmp <- NA
        dose_vals[i] <- tmp
      }
      dose_ins_row <- match(seq_along(labix), labix)
      curr_dose[dose_ins_row] <- as.numeric(dose_vals)
      values[[dat]][,'dose'] <- curr_dose
    })
  }
  factory_time <- function(dat, myid) {
    timeid <- sprintf('dt_%s_', myid)
    observeEvent({jj <- nrow(values[[dat]]); lapply(paste0(timeid, seq(jj)), function(i) input[[i]])}, {
      if(!values$init) return()
      labix <- values[[dat]][,'lab']
      curr_time <- unclass(values[[dat]][,'dt'])
      time_vals <- vector('list', length(labix))
      for(i in seq_along(time_vals)) {
        tmp <- input[[paste0(timeid, i)]]
        if(is.null(tmp)) tmp <- NA
        time_vals[[i]] <- tmp
      }
      time_ins_row <- match(seq_along(labix), labix)
      curr_time[time_ins_row] <- unlist(time_vals)
      curr_time <- as.POSIXct(curr_time, origin = '1970-01-01 00:00:00')
      values[[dat]][,'dt'] <- curr_time
    })
  }
  factory_ui <- function(dat, myid) {
    button1 <- sprintf('add_%s_button', myid)
#     button2 <- sprintf('delete_%s_button', myid)
    renderUI({
      if(!values$init) values$init <- TRUE
      ix <- isolate(c(TRUE, values[[dat]][,'valid'], TRUE))
      xx <- isolate(data2Inputs(values[[dat]], myid))
      # only add/delete should re-render
      input[[button1]]
#       input[[button2]]
      if(!is.null(xx)) do.call(tags$div, xx[ix])
    })
  }

  factory_add('infusion1dat', 'i1')
  factory_del('infusion1dat', 'i1')
  factory_dose('infusion1dat', 'i1')
  factory_time('infusion1dat', 'i1')
  factory_add('infusion2dat', 'i2')
  factory_del('infusion2dat', 'i2')
  factory_dose('infusion2dat', 'i2')
  factory_time('infusion2dat', 'i2')
  factory_add('infusion3dat', 'i3')
  factory_del('infusion3dat', 'i3')
  factory_dose('infusion3dat', 'i3')
  factory_time('infusion3dat', 'i3')
  factory_add('concdat', 'c')
  factory_del('concdat', 'c')
  factory_dose('concdat', 'c')
  factory_time('concdat', 'c')

  output$inf1 <- factory_ui('infusion1dat', 'i1')
  output$inf2 <- factory_ui('infusion2dat', 'i2')
  output$inf3 <- factory_ui('infusion3dat', 'i3')
  output$conc <- factory_ui('concdat', 'c')
  ####################    ####################    ####################    ####################    ####################

#   observe({
  observeEvent(input$runmodel, {
    def_start <- as.POSIXct("2022-01-01 00:00:00")
    def_inf <- data.frame(dose = 10, duration = 2, dt = def_start)
    def_con <- data.frame(dose = 10, dt = def_start + 16200) # 4.5 hours
    inf_duration <- 2
    in_infusmat1 <- checkInputMatrix(values$infusion1dat[values$infusion1dat[,'valid'],c("dose","dt")], def_inf, inf_duration)
    in_infusmat2 <- checkInputMatrix(values$infusion2dat[values$infusion2dat[,'valid'],c("dose","dt")], def_inf, inf_duration)
    in_infusmat3 <- checkInputMatrix(values$infusion3dat[values$infusion3dat[,'valid'],c("dose","dt")], def_inf, inf_duration)
    in_concmat <- checkInputMatrix(values$concdat[values$concdat[,'valid'],c("dose","dt")], def_con)

    #repair broken inputs
    usrwt <- as.numeric(input$wt)
    usralb <- as.numeric(input$alb)
    usrada <- as.numeric(input$ada)
    negNA <- function(x) is.na(x) || x < 0
    if(negNA(usrwt)) {
      usrwt <- 70
      updateTextInput(getDefaultReactiveDomain(), inputId = "wt", value = usrwt)
    }
    if(negNA(usralb)) {
      usrcreat <- 4
      updateTextInput(getDefaultReactiveDomain(), inputId = "alb", value = usralb)
    }
    if(negNA(usrada)) {
      usrada <- 0
      updateTextInput(getDefaultReactiveDomain(), inputId = "ada", value = usrada)
    }
    if(input$unit == 'mg') {
      # convert dose from "mg" to "mg/kg"
      in_infusmat1$dose <- in_infusmat1$dose / usrwt
      in_infusmat2$dose <- in_infusmat2$dose / usrwt
      in_infusmat3$dose <- in_infusmat3$dose / usrwt
    }

    ldat <- list(
      infList = list(in_infusmat1, in_infusmat2, in_infusmat3),
      conList = list(in_concmat)
    )

    myparams <- c(ldat, drug=input$drug, wt=usrwt, alb=usralb, ada=usrada)
    v$params <- myparams
    stuff <- do.call(setupModel, myparams)
    v$sched <- list(stuff[[1]][[1]], stuff[[1]][[2]], stuff[[1]][[3]])
    PKprof$Omega <- stuff[[5]]$Omega
    PKprof$p <- stuff[[5]]$p
    PKprof$wt <- stuff[[5]]$wt
    makePlots(stuff[[1]][[1]], stuff[[1]][[2]], stuff[[1]][[3]], stuff[[4]], stuff[[2]], stuff[[3]])
  })

  observeEvent(input$findTarget, {
    targetT <- as.numeric(input$targetTrough)
    if(!is.null(PKprof$p) && !is.na(targetT)) {
      tabdose <- as.numeric(input$ud)
      nDays <- 224 # number of days for 4 intervals of Q8WK, 7*8*4
      if(is.na(tabdose)) {
        # solve required dose
        dpd <- switch(input$ufrq,
                      Q1WK = 24 * 7 * 1,
                      Q2WK = 24 * 7 * 2,
                      Q3WK = 24 * 7 * 3,
                      Q4WK = 24 * 7 * 4,
                      Q5WK = 24 * 7 * 5,
                      Q6WK = 24 * 7 * 6,
                      Q7WK = 24 * 7 * 7,
                      Q8WK = 24 * 7 * 8,
                      24 * 7 * 2 # HERE
        )
        dopts <- seq(20)
        ptus <- vapply(dopts, function(i) getTrough(PKprof, i, dpd, nDays), numeric(2))
        dopt <- which.min(abs(targetT - ptus[2,]))
        updateTextInput(getDefaultReactiveDomain(), inputId = "ud", value = dopts[dopt])
      } else {
        # solve required frequency
        ptus <- vapply(1:8, function(i) getTrough(PKprof, tabdose, 24*7*i, nDays), numeric(2))
        dopt <- which.min(abs(targetT - ptus[2,]))
        dlab <- sprintf('Q%sWK', dopt)
        updateSelectInput(getDefaultReactiveDomain(), inputId = "ufrq", selected = dlab)
      }
    }
  })

  makePlots <- function(schedule1, schedule2, schedule3, concdat, tdrupper, tdrlower) {
    cols <- c("Individual" = "blue", 
              "Population" = "darkgoldenrod", 
              "Simulated" = "black",
              " Observed Drug Level" = "red")

    p1 <- bldplot(schedule1, 'Dosing 1', tdrupper, tdrlower) +
      geom_point(data=concdat, aes(x=yt/(24*7),y=y, color=colour), size = 3) +  # HERE rescale time here x=yt/(24*7) <- x=yt
      scale_colour_manual(
        values = cols,
        guide = guide_legend(
          override.aes = list(linetype = c("blank", rep("solid", 3)), shape = c(16, rep(NA, 3)))
        )
      ) +
      theme(plot.margin = ggplot2::margin(t=4,r=1,b=1,l=1, "lines")) +
      theme(legend.direction = "horizontal") +
      theme(legend.position = c(0.5, 1.2)) +
      labs(color= "")

    if(is.null(schedule2)) {
      p2 <- ggplot() + theme_void()
    } else {
      p2 <- bldplot(schedule2, 'Dosing 2', tdrupper, tdrlower) +
        scale_colour_manual(values = cols[1:3]) + theme(legend.position = "none")
    }
    if(is.null(schedule3)) {
      p3 <- ggplot() + theme_void()
    } else {
      p3 <- bldplot(schedule3, 'Dosing 3', tdrupper, tdrlower) +
        scale_colour_manual(values = cols[1:3]) + theme(legend.position = "none")
    }
    v$plot1 <- p1 + labs(title = 'dose schedule 1')
    v$plot2 <- p2 + labs(title = 'dose schedule 2')
    v$plot3 <- p3 + labs(title = 'dose schedule 3')
  }

  # user uploads CSV file
  observeEvent(v$dat, {
    if(!is.null(v$dat)) {
      ldat <- arrange_data(v$dat)
      # use `ldat` to update input form?
      usrwt <- as.numeric(ldat$covariates$wt)
      usralb <- as.numeric(ldat$covariates$alb)
      usrada <- as.numeric(ldat$covariates$ada)
      if(!is.null(usrwt)) updateTextInput(inputId='wt', value=ldat$covariates$wt)
      if(!is.null(usralb)) updateTextInput(inputId='alb', value=ldat$covariates$alb)
      if(!is.null(usrada)) updateTextInput(inputId='ada', value=ldat$covariates$ada)
      ldat$covariates <- NULL
      myparams <- c(ldat, drug=input$drug, wt=usrwt, alb=usralb, ada=usrada)
      v$params <- myparams
      stuff <- do.call(setupModel, myparams)
      v$sched <- list(stuff[[1]][[1]])
      PKprof$Omega <- stuff[[5]]$Omega
      PKprof$p <- stuff[[5]]$p
      PKprof$wt <- stuff[[5]]$wt
      makePlots(stuff[[1]][[1]], NULL, NULL, stuff[[4]], stuff[[2]], stuff[[3]])
    }
  })

  output$responseplot <- renderPlot(height = 600,{
    if(!is.null(v$plot1) && !is.null(v$plot2) && !is.null(v$plot3)) {
      grid.arrange(v$plot1, v$plot2, v$plot3, heights=c(5.5,4,4))
    } else {
      NULL
    }
  })

  tableInput <- function() {
    if(is.null(v$params)) return(NULL)

    i <- v$params[['infList']][[1]]
    c <- v$params[['conList']][[1]]
    w <- 1 # multiply by "weight" or not?
    c.s <- data.frame(Time = as.numeric(c[,'time']), Conc = c[,1], Amt = NA, Rate = NA, Duration = NA)
    # Amt = Rate * Duration -- is this correct?
    myrate <- as.numeric(i[,1]) / as.numeric(i[,3])
    i.s <- data.frame(Time = as.numeric(i[,'time']), Conc = NA, Amt = i[,1] * w, Rate = myrate, Duration = i[,3])
    t.out <- rbind(i.s,c.s)
    zero.amt <- which(!is.na(t.out[,'Amt']) & t.out[,'Amt']==0)
    if(length(zero.amt)) {
      t.out <- t.out[-zero.amt,]
    }
    t.out <- t.out[order(t.out[,'Time']),]
    t.out[,'Time'] <- (t.out[,'Time'] - min(t.out[,'Time']))
    t.out$weight <- v$params$wt
    t.out$alb <- v$params$alb
    t.out$ada <- v$params$ada
    t.out
  }

  output$profile1 <- DT::renderDT({
    tableInput()
  })

  output$peaktroughtable <- renderTable({
    n_ds <- length(v$sched)
    if(is.null(v$sched)) {
      pt_trgh <- matrix(NA, 3, 3)
      dose_plan <- c(5, 7.5, 10)
      n_ds <- 3
    } else {
      tp <- c(2, 6, 14) * 7 * 24
      sloops <- seq(min(3, n_ds))
      if(n_ds == 1) sloops <- c(1,1)
      pt_trgh <- t(vapply(v$sched[sloops], function(i) {
        i <- i[i$id == 'Individual',]
        i[i[,'t'] %in% tp, 'concentration']
      }, numeric(3)))
      dose_plan <- vapply(v$params[['infList']], function(i) i$dose[1], numeric(1))
    }
    colnames(pt_trgh) <- c('Week 2', 'Week 6', 'Week 14')
    # d <- sprintf('%s mg/kg', dose_plan)
    d <- paste0('Schedule ', seq(nrow(pt_trgh)))
    cbind(data.frame(Dosage = d), pt_trgh)[seq(n_ds),]
  })

  getTrough <- function(pk, d, f, nDays) {
    day10 <- nDays + 7
    inf.st.time <- seq(0, nDays*24, by = f)
    inf.e.time <- inf.st.time + 2
    tv <- seq(1, day10*24, by = 1) # HERE LC <- by = 48
    response20 <- getPatDat(clvec=pk$p[1],v1vec=pk$p[2],qvec=pk$p[3],v2vec=pk$p[4],
                            usrendt=day10*24, usrbdose=NULL, timevec=tv,
                            usrbt=NULL, usrwt=pk$wt,
                            usriet=inf.e.time, usrist=inf.st.time,
                            usridose=rep(d, length(inf.st.time)), Omega=pk$Omega)
    ss_start <- 24 * 7 * 8 * 3 # HERE LC   24 * 7 * 8 * 3: number of hours for 3 intervals of Q8WK
    getConcRange(response20, (ss_start + 1), (ss_start + f) ) # HERE LC: 1 dosing time interval after 3 of Q8WK (i.e., 168 days after)
  }

  output$utable <- renderTable({
    if(is.null(PKprof$p)) {
      ptu <- rep(NA, 2)
    } else {
      tabdose <- as.numeric(input$ud)
      targetT <- as.numeric(input$targetTrough)
      if(!is.na(tabdose) && tabdose < 0) {
        tabdose <- 10
        updateTextInput(getDefaultReactiveDomain(), inputId = "ud", value=tabdose)
      }
      dpd <- switch(input$ufrq,
                    Q1WK = 24 * 7 * 1,
                    Q2WK = 24 * 7 * 2,
                    Q3WK = 24 * 7 * 3,
                    Q4WK = 24 * 7 * 4,
                    Q5WK = 24 * 7 * 5,
                    Q6WK = 24 * 7 * 6,
                    Q7WK = 24 * 7 * 7,
                    Q8WK = 24 * 7 * 8,
                    24 * 7 * 2 # HERE
      )
      nDays <- 224 # HERE LC number of days for 4 intervals of Q8WK
      ptu <- getTrough(PKprof, tabdose, dpd, nDays)
    }
    d <- c("Trough (mcg/mL)") #  HERE LC from d <- c("Peak (mcg/mL)", "Trough (mcg/mL)")
    udf <- data.frame(Col1 = d, Col2 = round(ptu[2], 2)) #  HERE LC from Col2 = round(ptu, 2)
    names(udf)<- c(" "," ")
    udf
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { sprintf('TDM_plot_%s.pdf', Sys.Date()) },
    content = function(file, width = 8, height = 8) {
      pdf(file, width = width, height = height, onefile = TRUE)
      grid.arrange(v$plot1, v$plot2, v$plot3, heights = c(5, 4, 4))
      dev.off()
  })

  output$downloadTable <- downloadHandler(
    filename = function() { sprintf('TDM_profile_%s.csv', Sys.Date()) },
    content = function(file) {
      write.csv(tableInput(), file, row.names = FALSE)
  })

  output$downloadTable1 <- downloadHandler(
    filename = function() { sprintf('TDM_profile_%s.csv', Sys.Date()) },
    content = function(file) {
      write.csv(tableInput(), file, row.names = FALSE)
  })

  output$imgdat <- renderUI({
    tags$img(src = "https://raw.githubusercontent.com/michaelleewilliams/michaelleewilliams.github.io/master/pictures/datex.png",height="50%", width="50%", align="center")
  })

  output$imgscn <- renderUI({
    tags$img(src = "https://raw.githubusercontent.com/michaelleewilliams/michaelleewilliams.github.io/master/pictures/scn.png",height="100%", width="100%", align="center")
  })
}

shinyApp(ui, server)
