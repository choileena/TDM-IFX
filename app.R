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
library(deSolve)
library(cowplot)
library(plotly)
library(gapminder)
library(shinyTime)
library(tidyverse)

# Load each helpers file into its own environment so internal functions
# (getpatdat, simulate_*, pkprof_est_*, etc.) don't overwrite each other

env_mtx <- new.env(parent = globalenv())
sys.source('helpers_temp_mtx_v9.R',   envir = env_mtx)

env_hui <- new.env(parent = globalenv())
sys.source('helpers_temp_hui_v8_5.R', envir = env_hui)

# Wrap each setupModel so its entire internal call chain runs in the right environment
setupModel_mtx <- function(...) {
  args <- list(...)
  environment(env_mtx$setupModel) <- env_mtx
  do.call(env_mtx$setupModel, args)
}
setupModel_hui <- function(...) {
  args <- list(...)
  environment(env_hui$setupModel) <- env_hui
  do.call(env_hui$setupModel, args)
}

makePlots_mtx <- env_mtx$makePlots
makePlots_hui <- env_hui$makePlots

# ── data2Inputs (shared) ──────────────────────────────────────────────────────
data2Inputs <- function(dd, myid) {
  doseLabel <- c('Bolus', 'Infusion', 'Conc')[match(substr(myid, 1, 1), c('b', 'i', 'c'))]
  nn        <- nrow(dd)
  adder     <- sprintf('Shiny.onInputChange( \"add_%s_button\" , this.id, {priority: \"event\"})', myid)
  add_label <- if (doseLabel == 'Conc') 'add conc' else 'add dose'
  addit     <- actionButton('addit', label = add_label, icon = icon('plus'), onclick = adder,
                            style = 'padding:4px;font-size:80%;margin-top:4px;')
  
  header_cells <- tagList(
    tags$th(doseLabel,     style = "width:65px;"),
    tags$th("Date & Time", style = "width:190px;")
  )
  if (doseLabel != 'Conc') {
    header_cells <- tagList(header_cells, tags$th("Dur (Hrs)", style = "width:70px;"))
  }
  header_row <- tags$tr(header_cells)
  
  data_rows <- lapply(seq_len(nn), function(i) {
    dt_value   <- as.POSIXct(dd[i, 'dt'], origin = "1970-01-01")
    dose_input <- textInput(paste0('dose_', myid, '_', i), label = NULL,
                            value = as.character(dd[i, 'dose']), width = '60px')
    date_input <- airDatepickerInput(paste0('date_', myid, '_', i), label = NULL,
                                     timepicker = TRUE, dateFormat = "yyyy-MM-dd",
                                     timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                     value = dt_value, width = '180px', readonly = TRUE)
    cells <- tagList(tags$td(dose_input), tags$td(date_input))
    if (doseLabel != 'Conc') {
      dur_input <- textInput(paste0('duration_', myid, '_', i), label = NULL,
                             value = as.character(dd[i, 'duration']), width = '60px')
      cells <- tagList(cells, tags$td(dur_input))
    }
    tags$tr(cells)
  })
  
  tags$div(
    tags$table(class = "input-table", header_row, tagList(data_rows)),
    tags$div(addit)
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  useShinyjs(),
  
  tags$head(tags$style(HTML("
    .disclaimer-banner {
      background-color: #fff5f5; color: #8B0000;
      border-left: 6px solid #8B0000; padding: 14px 18px;
      margin-bottom: 20px; border-radius: 6px;
      font-size: 15px; line-height: 1.45;
    }
    .input-table { border-collapse: separate; border-spacing: 4px 0; }
    .input-table th { font-weight: normal; font-size: 13px; padding-bottom: 2px; white-space: nowrap; }
    .input-table td { vertical-align: top; padding-right: 4px; }
    .app-selector { padding: 10px 0 5px 0; }
  "))),
  
  tags$script(HTML("
    $(document).on('change', '.air-datepicker-input', function() {
      $(this).trigger('change');
    });
  ")),
  
  tags$div(
    class = "disclaimer-banner",
    tags$strong("Clinical Disclaimer: "),
    "This application is intended for research and informational purposes only. ",
    "Predictions and estimates may not accurately reflect individual patient circumstances, local patient populations, or institutional practices. ",
    "Model performance may vary across clinical settings and patient populations. ",
    "This tool is not intended to replace professional clinical judgment, institutional protocols, or medical advice. ",
    "All clinical and dosing decisions remain the responsibility of the treating clinician. ",
    "Use of this application is at the user's own risk."
  ),
  
  titlePanel("Therapeutic Drug Monitoring"),
  
  # ── App selector dropdown ──
  fluidRow(
    column(3,
           tags$div(class = "app-selector",
                    selectInput("active_app", label = strong("Select Model:"),
                                choices = c("Taylor App (Blackman & Taylor)" = "mtx",
                                            "Hui App (Blackman & Hui)"   = "hui"),
                                width = "280px")
           )
    )
  ),
  
  tabsetPanel(
    tabPanel("Application",
             fluidRow(
               column(3,
                      hr(),
                      fileInput('file1', 'Upload Dosing Profile (CSV)'),
                      hr(),
                      helpText(strong("Drug: "), "High-Dose Methotrexate"),
                      hr(),
                      tabsetPanel(type = "pills",
                                  tabPanel("Baseline Information",
                                           helpText(em("Input baseline patient characteristics.")),
                                           textInput("bsa", "BSA (m²)", value = "1.97"),
                                           selectInput("sex", "Sex:", c("Male", "Female")),
                                           # Hui-only fields (hidden for MTX app)
                                           conditionalPanel("input.active_app == 'hui'",
                                                            textInput("height_cm", "Height (cm)", value = "170"),
                                                            airDatepickerInput("dob", label = "Date of Birth",
                                                                               value = Sys.Date() - 365 * 10,
                                                                               maxDate = Sys.Date(),
                                                                               dateFormat = "yyyy-MM-dd",
                                                                               width = "180px", readonly = TRUE)
                                           ),
                                           fluidRow(
                                             column(12, textInput("SCR_mgdl1", "Baseline Serum Creatinine (mg/dL)", value = "0.77")),
                                             column(12, shinyjs::hidden(airDatepickerInput("SCR_time1", label = "Time",
                                                                                           timepicker = TRUE, readonly = TRUE, dateFormat = "yyyy-MM-dd",
                                                                                           timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                                                                           value = Sys.time())))
                                           )
                                  ),
                                  tabPanel("Dosing",
                                           helpText(em("Provide total amount administered per infusion in mg as well as duration in hours. Leave amount blank to cancel an extra dose.")),
                                           htmlOutput('inf1')
                                  ),
                                  tabPanel("Drug Levels",
                                           helpText(em("To generate population-level predictions only, leave the drug level amount blank and click 'Create Output'.")),
                                           htmlOutput('conc')
                                  ),
                                  tabPanel("Serum Creatinine Levels Post-infusion",
                                           helpText(em("Add up to 5 additional measured serum creatinine levels. Be sure to set the time of each observed level.")),
                                           fluidRow(
                                             column(12, textInput("SCR_mgdl2", "Serum Creatinine 2 (mg/dL)", value = "0.80")),
                                             column(12, airDatepickerInput("SCR_date2", label = "Date & Time 2", timepicker = TRUE,
                                                                           readonly = TRUE, dateFormat = "yyyy-MM-dd",
                                                                           timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                                                           value = Sys.time() + 86400, width = "180px"))
                                           ),
                                           fluidRow(
                                             column(12, textInput("SCR_mgdl3", "Serum Creatinine 3 (mg/dL)", value = "")),
                                             column(12, airDatepickerInput("SCR_date3", label = "Date & Time 3", timepicker = TRUE,
                                                                           readonly = TRUE, dateFormat = "yyyy-MM-dd",
                                                                           timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                                                           value = Sys.time() + 86400, width = "180px"))
                                           ),
                                           fluidRow(
                                             column(12, textInput("SCR_mgdl4", "Serum Creatinine 4 (mg/dL)", value = "")),
                                             column(12, airDatepickerInput("SCR_date4", label = "Date & Time 4", timepicker = TRUE,
                                                                           readonly = TRUE, dateFormat = "yyyy-MM-dd",
                                                                           timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                                                           value = Sys.time() + 86400, width = "180px"))
                                           ),
                                           fluidRow(
                                             column(12, textInput("SCR_mgdl5", "Serum Creatinine 5 (mg/dL)", value = "")),
                                             column(12, airDatepickerInput("SCR_date5", label = "Date & Time 5", timepicker = TRUE,
                                                                           readonly = TRUE, dateFormat = "yyyy-MM-dd",
                                                                           timepickerOpts = timepickerOptions(timeFormat = "HH:mm"),
                                                                           value = Sys.time() + 86400, width = "180px"))
                                           )
                                  )
                      )
               ),
               column(6,
                      tags$br(),
                      actionButton('runmodel', 'Create Output'),
                      plotlyOutput("responseplot", width = "60vw", height = "70vh"),
                      align = 'center',
                      helpText(strong("Recommendation:")),
                      conditionalPanel("input.active_app == 'mtx'",
                                       helpText("Follow population prediction line until first observed drug level then follow individual Taylor model prediction line")
                      ),
                      conditionalPanel("input.active_app == 'hui'",
                                       helpText("Follow population prediction line until first observed drug level then follow individual Hui model prediction line")
                      ),
                      hr()
               )
             )
    ),
    tabPanel("Documentation",
             h1("Functionality"),
             p("This application takes as input patient characteristics/dosing schedule and presents the patient's predicted concentration-time curve. The concentration predictions are
             based on the solutions to the standard three-compartment pharmacometric (PK) ordinary differential equations. The predictions include population-level estimates and 
        will also include individual-level expected response when an observed blood level draw is available. The individual-level estimates are rendered with empirical Bayesian estimates (EBEs),
        so they require the observed blood level draw for estimation. The Bayesian process is based on the minimization of the likelihood with respect to the individual
        random effects [1]. The Taylor App uses Blackman et al. [2] and Taylor et al. [3] models; the Hui App uses Blackman et al. and Hui et al. [4] models."),
             h1("Resources"),
             HTML("<p>The models implemented in this application are all available in the scientific literature and cited below.
           The calculation of the EBEs for individual random effects is outlined in
           <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/'>Kang 2012</a> [1]
           and the code to calculate the EBEs is adapted from the open source TDM software
           <a href='https://mrgsolve.org/'>mrgsolve</a>.
           </p>"),
             h1("Bibliography"),
             tags$ol(
               tags$li("Kang, D., Bae, K.-S., Houk, B. E., Savic, R. M. & Karlsson, M. O. Standard Error of Empirical Bayes Estimate in NONMEM® VI. Korean J Physiol Pharmacol (2012)."),
               tags$li("Blackman et al. Development and Validation of High-Dose Methotrexate Population Pharmacokinetic Models to Inform Clinical Decisions on Dosing, European Journal of Clinical Pharmacology (2026)."),
               tags$li("Taylor et al. MTXPK.org: A Clinical Decision Support Tool Evaluating High-Dose Methotrexate Pharmacokinetics to Inform Post-Infusion Care and Use of Glucarpidase, Clinical Pharmacology and Pharmacometrics (2020)."),
               tags$li("Hui et al. Population Pharmacokinetic Study and Individual Dose Adjustments of High-Dose Methotrexate in Chinese Pediatric Patients With Acute Lymphoblastic Leukemia or Osteosarcoma, The Journal of Clinical Pharmacology (2018).")
             )
    ),
    tabPanel("View Profile",
             DT::dataTableOutput('profile1'),
             downloadButton('downloadTable1', 'Download Data (CSV)')
    ),
    tabPanel("Predictions",
             DT::dataTableOutput('profile2'),
             downloadButton('downloadTable2', 'Download Data (CSV)')
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  showModal(modalDialog(
    title = "Clinical Disclaimer",
    div(style = "color:#8B0000; line-height:1.5;",
        strong("This tool provides estimates only. "),
        "Predictions may not accurately reflect your patient population or individual patient circumstances. ",
        "This tool is not a substitute for independent clinical judgment, institutional protocols, or applicable standards of care. ",
        "All dosing decisions remain the responsibility of the treating clinician. Use at your own risk."
    ),
    easyClose = TRUE, footer = modalButton("I Understand"), size = "m"
  ))
  
  v      <- reactiveValues(dat = NULL, params = NULL, sched = NULL)
  PKprof <- reactiveValues(p = NULL, Omega = NULL, unit = 'mg')
  
  observeEvent(input$file1, {
    if (!is.null(input$file1)) {
      inFile <- input$file1
      if (!is.null(inFile) && inFile$type %in% c('text/csv', 'text/comma-separated-values', 'text/plain'))
        v$dat <- read.csv(inFile$datapath)
      else
        v$dat <- NULL
    }
  })
  
  def_time1 <- Sys.time()
  
  values <- reactiveValues(
    infusion1dat = data.frame(dose = 2330, dt = def_time1,         lab = 1, duration = 24, valid = TRUE),
    concdat      = data.frame(dose = 20,   dt = def_time1 + 86400, lab = 1,               valid = TRUE),
    init = FALSE
  )
  
  # ── Factories ──────────────────────────────────────────────────────────────
  factory_add <- function(dat, myid) {
    button <- sprintf('add_%s_button', myid)
    observeEvent(input[[button]], {
      if (!values$init) return()
      ix <- which(values[[dat]][, 'valid'])
      if (length(ix) == 0) {
        next_time <- as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:00"))
        new_duration <- 24
      } else {
        first_time   <- as.POSIXct(values[[dat]][1, 'dt'])
        first_dur    <- as.numeric(values[[dat]][1, 'duration'])
        next_time    <- first_time + first_dur * 3600
        new_duration <- max(0, 24 - first_dur)
      }
      values[[dat]] <- rbind(values[[dat]],
                             data.frame(dose = NA, dt = next_time, lab = nrow(values[[dat]]) + 1,
                                        duration = new_duration, valid = TRUE))
    })
  }
  
  factory_addconc <- function(dat, myid) {
    button <- sprintf('add_%s_button', myid)
    observeEvent(input[[button]], {
      if (!values$init) return()
      ix <- which(values[[dat]][, 'valid'])
      next_time <- if (length(ix) == 0)
        as.POSIXct(format(Sys.time(), "%Y-%m-%d %H:00"))
      else
        as.POSIXct(format(max(values[[dat]][ix, 'dt']), "%Y-%m-%d %H:00")) + 3600 * 12
      values[[dat]] <- rbind(values[[dat]],
                             data.frame(dose = NA, dt = next_time, lab = nrow(values[[dat]]) + 1, valid = TRUE))
    })
  }
  
  factory_dose <- function(dat, myid) {
    doseid <- sprintf('dose_%s_', myid)
    observeEvent({
      jj <- nrow(values[[dat]])
      lapply(paste0(doseid, seq(jj)), function(i) input[[i]])
    }, {
      if (!values$init) return()
      labix     <- values[[dat]][, 'lab']
      curr_dose <- values[[dat]][, 'dose']
      dose_vals <- character(length(labix))
      for (i in seq_along(dose_vals)) {
        tmp <- input[[paste0(doseid, i)]]
        if (is.null(tmp)) tmp <- NA
        dose_vals[i] <- tmp
      }
      curr_dose[match(seq_along(labix), labix)] <- as.numeric(dose_vals)
      values[[dat]][, 'dose'] <- curr_dose
    })
  }
  
  factory_time <- function(dat, myid) {
    observe({
      if (!values$init) return()
      labix   <- values[[dat]][, 'lab']
      curr_dt <- values[[dat]][, 'dt']
      changed <- FALSE
      for (i in seq_along(labix)) {
        dt_val <- input[[paste0('date_', myid, '_', i)]]
        if (is.null(dt_val)) next
        new_dt <- tryCatch(as.POSIXct(dt_val), error = function(e) NULL)
        if (is.null(new_dt) || is.na(new_dt)) next
        row <- match(i, labix)
        if (!is.na(row) && !identical(curr_dt[row], new_dt)) {
          curr_dt[row] <- new_dt
          changed <- TRUE
        }
      }
      if (changed) values[[dat]][, 'dt'] <- curr_dt
    })
  }
  
  factory_duration <- function(dat, myid) {
    durationid <- sprintf('duration_%s_', myid)
    observeEvent({
      jj <- nrow(values[[dat]])
      lapply(paste0(durationid, seq(jj)), function(i) input[[i]])
    }, {
      if (!values$init) return()
      labix    <- values[[dat]][, 'lab']
      curr_dur <- values[[dat]][, 'duration']
      dur_vals <- character(length(labix))
      for (i in seq_along(dur_vals)) {
        tmp <- input[[paste0(durationid, i)]]
        if (is.null(tmp)) tmp <- NA
        dur_vals[i] <- tmp
      }
      curr_dur[match(seq_along(labix), labix)] <- as.numeric(dur_vals)
      values[[dat]][, 'duration'] <- curr_dur
    })
  }
  
  factory_ui <- function(dat, myid) {
    button1 <- sprintf('add_%s_button', myid)
    renderUI({
      if (!values$init) values$init <- TRUE
      input[[button1]]
      isolate(data2Inputs(values[[dat]][values[[dat]][, 'valid'], ], myid))
    })
  }
  
  factory_add('infusion1dat', 'i1')
  factory_dose('infusion1dat', 'i1')
  factory_time('infusion1dat', 'i1')
  factory_duration('infusion1dat', 'i1')
  
  factory_addconc('concdat', 'c')
  factory_dose('concdat', 'c')
  factory_time('concdat', 'c')
  
  output$inf1 <- factory_ui('infusion1dat', 'i1')
  output$conc <- factory_ui('concdat', 'c')
  
  # ── Shared SCR/input helpers ───────────────────────────────────────────────
  get_scr_vals <- function() {
    sapply(1:5, function(i) {
      x <- input[[paste0("SCR_mgdl", i)]]
      if (is.null(x) || x == "") NA_real_ else as.numeric(x)
    })
  }
  
  get_scr_times <- function() {
    c(
      {
        dose_dt <- values$infusion1dat[values$infusion1dat[, 'valid'], 'dt']
        if (length(dose_dt) == 0 || is.na(dose_dt[1])) as.numeric(Sys.time())
        else as.numeric(as.POSIXct(dose_dt[1]))
      },
      sapply(2:5, function(i) {
        dt_val <- input[[paste0("SCR_date", i)]]
        if (is.null(dt_val)) return(NA_real_)
        as.numeric(as.POSIXct(dt_val))
      })
    )
  }
  
  get_common_inputs <- function(env) {
    def_start    <- Sys.time()
    def_inf      <- data.frame(dose = 2330, duration = 24, dt = def_start)
    def_con      <- data.frame(dose = 20,   dt = def_start + 86400)
    in_infusmat1 <- env$checkInputMatrix(
      values$infusion1dat[values$infusion1dat[, 'valid'], c("dose", "duration", "dt")],
      def_inf, 24)
    raw_concdat  <- values$concdat[values$concdat[, 'valid'], ]
    has_conc     <- nrow(raw_concdat) > 0 &&
      any(!is.na(suppressWarnings(as.numeric(raw_concdat$dose))) &
            suppressWarnings(as.numeric(raw_concdat$dose)) > 0)
    in_concmat   <- env$checkInputMatrix(raw_concdat[, c("dose", "dt")], def_con)
    bsa          <- as.numeric(input$bsa)
    if (is.na(bsa) || bsa < 0) { bsa <- 1.97; updateTextInput(session, "bsa", value = bsa) }
    pt_gender    <- input$sex
    scr_vals     <- get_scr_vals()
    if (is.na(scr_vals[1]) || scr_vals[1] < 0) {
      scr_vals[1] <- 68.08 / 88.4
      updateTextInput(session, "SCR_mgdl1", value = scr_vals[1])
    }
    scr_times <- get_scr_times()
    list(in_infusmat1 = in_infusmat1, in_concmat = in_concmat, has_conc = has_conc,
         bsa = bsa, pt_gender = pt_gender, scr_vals = scr_vals, scr_times = scr_times)
  }
  
  # ── Plotly post-processing (shared) ───────────────────────────────────────
  style_plotly <- function(gp, app_type) {
    second_model_label <- if (app_type == "mtx") "Taylor model" else "Hui model"
    second_model_regex <- if (app_type == "mtx") "Taylor model"  else "Hui model"
    
    for (i in seq_along(gp$x$data)) {
      if (!is.null(gp$x$data[[i]]$name) && grepl("Serum Creatinine", gp$x$data[[i]]$name))
        gp$x$data[[i]]$yaxis <- "y2"
    }
    
    for (i in seq_along(gp$x$data)) {
      raw_name <- gp$x$data[[i]]$name
      if (is.null(raw_name)) { gp$x$data[[i]]$showlegend <- FALSE; next }
      clean_name <- gsub("^\\(|\\,1\\)$", "", raw_name)
      clean_name <- gsub("Blackman et al\\.", "Blackman model", clean_name)
      clean_name <- gsub("MTXPK\\.org",      "Taylor model",   clean_name)
      clean_name <- gsub("Hui et al\\.",      "Hui model",      clean_name)
      
      if (clean_name == " ") {
        gp$x$data[[i]]$name           <- " "
        gp$x$data[[i]]$showlegend     <- TRUE
        gp$x$data[[i]]$visible        <- TRUE
        gp$x$data[[i]]$marker$opacity <- 0
        gp$x$data[[i]]$hoverinfo      <- "none"
        next
      }
      if (clean_name == "" || grepl("Trace", clean_name) || grepl("df_ribbon", clean_name)) {
        gp$x$data[[i]]$showlegend <- FALSE
      } else {
        gp$x$data[[i]]$name       <- clean_name
        gp$x$data[[i]]$showlegend <- TRUE
        is_point <- clean_name %in% c("Observed Drug Level", "Serum Creatinine (mg/dL)", "Glucarpidase Consensus Guidelines")
        if (is_point) {
          gp$x$data[[i]]$mode       <- "markers"
          gp$x$data[[i]]$line$width <- 0
        } else {
          gp$x$data[[i]]$mode <- "lines"
          if (grepl(second_model_regex, clean_name)) {
            gp$x$data[[i]]$line$dash  <- "dash"
            gp$x$data[[i]]$line$width <- 1.5
          } else {
            gp$x$data[[i]]$line$dash <- "solid"
          }
          if (grepl("Blackman model", clean_name))
            gp$x$data[[i]]$line$width <- if (!grepl("probability", clean_name)) 2.5 else 1
          else
            gp$x$data[[i]]$line$width <- 1.5
        }
        if (grepl("Serum Creatinine", clean_name)) gp$x$data[[i]]$yaxis <- "y2"
        gp$x$data[[i]]$legendgroup <- clean_name
      }
    }
    gp
  }
  
  apply_layout <- function(gp, tick_positions, plot_range, level_order) {
    trace_names   <- sapply(gp$x$data, function(x) x$name)
    order_indices <- match(level_order, trace_names)
    order_indices <- order_indices[!is.na(order_indices)]
    other_indices <- setdiff(seq_along(gp$x$data), order_indices)
    gp$x$data     <- gp$x$data[c(order_indices, other_indices)]
    
    gp %>% layout(
      legend  = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2, traceorder = "normal"),
      margin  = list(r = 80),
      yaxis   = list(range = plot_range, tickmode = "auto", ticks = "outside",
                     ticklen = 5, tickwidth = 1, tickcolor = "black"),
      yaxis2  = list(overlaying = "y", side = "right", range = plot_range,
                     tickvals = tick_positions, ticktext = round(tick_positions / 30, 2),
                     title = list(text = "Serum Creatinine (mg/L)", font = list(size = 14), color = "purple"),
                     tickfont = list(size = 12, color = "purple"),
                     ticks = "outside", ticklen = 5, tickwidth = 1, tickcolor = "black")
    )
  }
  
  # ── Run model ──────────────────────────────────────────────────────────────
  observeEvent(input$runmodel, {
    showModal(modalDialog(title = NULL, "Calculating, please wait...",
                          footer = NULL, easyClose = FALSE))
    
    #ci <- get_common_inputs()
    app_type <- input$active_app
    ci <- get_common_inputs(if (app_type == "mtx") env_mtx else env_hui)
    if (app_type == "mtx") {
      myparams <- c(ci, list(drug = input$drug))
      stuff    <- do.call(setupModel_mtx, myparams)
    } else {
      height_cm <- as.numeric(input$height_cm)
      if (is.na(height_cm) || height_cm < 0) { height_cm <- 170; updateTextInput(session, "height_cm", value = 170) }
      myparams <- c(ci, list(drug = input$drug, height_cm = height_cm, dob = input$dob))
      stuff    <- do.call(setupModel_hui, myparams)
    }
    
    v$params <- c(myparams, list(app_type = app_type))
    v$sched  <- list(stuff[[1]], stuff[[2]])
    PKprof$Omega <- stuff[[3]]$Omega
    PKprof$p     <- stuff[[3]]
    
    removeModal()
    
    output$responseplot <- renderPlotly({
      df      <- stuff[[2]]
      concdat <- stuff[[1]]
      make_fn <- if (app_type == "mtx") makePlots_mtx else makePlots_hui
      p       <- make_fn(stuff[[2]], stuff[[1]])
      gp      <- ggplotly(p)
      
      con1  <- as.numeric(unlist(df[, 'Conc']))
      con2  <- c(50, 30, 10, 5) * 0.454
      con3  <- if (nrow(concdat) > 0) concdat[, 'y'] else numeric(0)
      con4  <- as.numeric(unlist(df[, 'SCR_mmol'])) / 88.4
      myvec <- c(con1, con2, unlist(con3), con4)
      myvec <- myvec[is.finite(myvec)]
      y_range        <- c(-0.5, max(myvec) + 2)
      ceiling_val    <- 5 * round(y_range[2] / 5)
      tick_positions <- if (ceiling_val >= 30) seq(0, ceiling_val, by = 10) else seq(0, ceiling_val, by = 5)
      plot_range     <- c(-0.5, y_range[2] * 1.05)
      
      gp <- style_plotly(gp, app_type)
      
      if (app_type == "mtx") {
        level_order <- c("Observed Drug Level", " ", "Serum Creatinine (mg/L)",
                         "Glucarpidase Consensus Guidelines",
                         "Population level: Blackman model", "Individual level: Blackman model",
                         "Population level: Taylor model",   "Individual level: Taylor model",
                         "Upper 95% population probability: Blackman model",
                         "Lower 95% population probability: Blackman model")
      } else {
        level_order <- c("Observed Drug Level", " ", "Serum Creatinine (mg/L)",
                         "Glucarpidase Consensus Guidelines",
                         "Population level: Blackman model", "Individual level: Blackman model",
                         "Population level: Hui model",      "Individual level: Hui model",
                         "Upper 95% population probability: Blackman model",
                         "Lower 95% population probability: Blackman model")
      }
      
      apply_layout(gp, tick_positions, plot_range, level_order)
    })
  })
  
  # ── Tables ─────────────────────────────────────────────────────────────────
  tableInput <- function() {
    if (is.null(v$params)) return(NULL)
    i      <- v$params[["in_infusmat1"]]
    c      <- v$params[["in_concmat"]]
    myrate <- round(as.numeric(i[, 1]) / as.numeric(i[, "duration"]), 2)
    i.s    <- data.frame(Time = as.numeric(i[, 'time']) / 3600, DV = NA, AMT = i[, 1],
                         MDV = 1, EVID = 1, Rate = myrate, Duration = i[, "duration"], SCR = NA)
    if (isTRUE(v$params[["has_conc"]])) {
      c.s   <- data.frame(Time = as.numeric(c[, 'time']) / 3600, DV = c[, 'dose'], AMT = NA,
                          MDV = 0, EVID = 0, Rate = NA, Duration = NA, SCR = NA)
      t.out <- rbind(i.s, c.s)
    } else {
      t.out <- i.s
    }
    t.out <- t.out[!(!is.na(t.out$AMT) & t.out$AMT == 0), ]
    t.out <- t.out[order(t.out$Time), ]
    
    scr_vals  <- round(v$params$scr_vals, 2)
    scr_times <- v$params$scr_times / 3600
    t.out[which(t.out$Time == min(t.out$Time)), "SCR"] <- scr_vals[1]
    
    if (length(scr_times) > 1) {
      scr.s  <- data.frame(Time = scr_times[2:length(scr_times)], DV = NA, AMT = NA,
                           MDV = 0, EVID = 0, Rate = NA, Duration = NA,
                           SCR = scr_vals[2:length(scr_vals)])
      t.out1 <- rbind(t.out, scr.s)
    } else {
      t.out1 <- t.out
    }
    t.out1 <- t.out1[order(t.out1$Time), ]
    t.out1$Time <- round(t.out1$Time - min(t.out$Time), 2)
    
    dups <- as.numeric(names(which(table(t.out1$Time) > 1)))
    for (tt in dups) {
      idx <- which(t.out1$Time == tt)
      t.out1[idx[1], "DV"]  <- t.out1[idx, "DV"][!is.na(t.out1[idx, "DV"])][1]
      t.out1[idx[1], "SCR"] <- t.out1[idx, "SCR"][!is.na(t.out1[idx, "SCR"])][1]
      t.out1 <- t.out1[-idx[-1], ]
    }
    t.out1$bsa <- v$params$bsa
    t.out1$Sex <- v$params$pt_gender
    t.out1 %>%
      filter(!is.na(Time)) %>%
      rename(`BSA (m²)` = bsa, `SCR (mg/dL)` = SCR) %>%
      fill(`SCR (mg/dL)`, .direction = "down")
  }
  
  table2Input <- function() {
    if (is.null(v$sched[[2]])) return(NULL)
    d   <- v$sched[[2]]
    dat <- d %>%
      filter(time %in% seq(0, 600, by = 6)) %>%
      dplyr::select(time, bsa, SCR_mmol, pt_gender, Conc, tag) %>%
      mutate(pt_gender = ifelse(pt_gender %in% c(0, "0"), "Male",
                                ifelse(pt_gender %in% c(1, "1"), "Female", as.character(pt_gender))),
             Conc     = round(Conc, 2),
             SCR_mmol = round(SCR_mmol / 88.4, 2)) %>%
      rename(`BSA (m²)` = bsa, `SCR (mg/dL)` = SCR_mmol,
             `Concentration (mg/L)` = Conc, Sex = pt_gender, `Time (hrs)` = time)
    dat %>% pivot_wider(
      id_cols    = c(`Time (hrs)`, `BSA (m²)`, `SCR (mg/dL)`, Sex),
      names_from = tag, values_from = `Concentration (mg/L)`
    ) %>% 
    relocate(`Lower 95% population probability: Blackman et al.`,
             .before = `Upper 95% population probability: Blackman et al.`)
  }
  
  output$profile1 <- DT::renderDT({ tableInput() })
  output$profile2 <- DT::renderDT({ table2Input() })
  
  output$downloadTable1 <- downloadHandler(
    filename = function() sprintf('TDM_profile_%s.csv',     Sys.Date()),
    content  = function(file) write.csv(tableInput(),  file, row.names = FALSE))
  output$downloadTable2 <- downloadHandler(
    filename = function() sprintf('TDM_predictions_%s.csv', Sys.Date()),
    content  = function(file) write.csv(table2Input(), file, row.names = FALSE))
}

shinyApp(ui, server)