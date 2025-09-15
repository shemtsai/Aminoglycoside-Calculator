# Package management
required_packages <- c("shiny", "ggplot2", "shinydashboard", "shinyTime", "bslib")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

library(shiny)
library(ggplot2)
library(shinydashboard)
library(shinyTime)
library(bslib)

# Theme
bs_theme <- bs_theme(
    version = 5,
    bootswatch = "minty",
    primary = "#2563eb",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    "border-radius" = "0.75rem"
)

ui <- fluidPage(
    theme = bs_theme,
    tags$head(
        tags$style(HTML("
            .card { border-radius: 18px !important; box-shadow: 0 4px 24px rgba(37,99,235,0.08); }
            .btn-primary { background-color: #2563eb; border-color: #2563eb; font-weight: 500; border-radius: 8px; }
            .form-control, .selectize-input { border-radius: 8px; }
            label, .control-label { font-weight: 600; }
            .nav-pills .nav-link.active { background: #2563eb; }
        "))
    ),
    titlePanel(
        div(
            icon("flask"),
            span("Aminoglycoside Calculator", style = "font-weight:600; letter-spacing:1px;")
        )
    ),
    tabsetPanel(
        tabPanel("Standard PK Calculator",
            fluidRow(
                column(4,
                    wellPanel(
                        radioButtons("input_method", "Select Input Method:",
                            choices = c("Enter Date/Time and Dose" = "date_time_dose", 
                            "Enter Concentrations and Times" = "concentration_time"
                                        ),
                            inline = TRUE
                        ),

                        conditionalPanel(
                            condition = "input.input_method == 'concentration_time'",
                            numericInput("c1", "First Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
                            numericInput("c2", "Second Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1),
                            numericInput("t1", "Time of First Level (t₁, hours) after START of infusion:", value = 4, min = 0, step = 0.1),
                            numericInput("t2", "Time of Second Level (t₂, hours) after START of infusion:", value = 12, min = 0, step = 0.1)
                        ),
                        conditionalPanel(
                            condition = "input.input_method == 'date_time_dose'",
                            dateInput("dose_date", "Date of Dose Administration (mm:dd):", value = Sys.Date(), format = "mm/dd"),
                            timeInput("dose_time", "Time of Dose Administration (hh:mm):", value = strptime("08:00", format="%H:%M"), seconds = FALSE),
                            dateInput("date1", "Date of First Level (t₁):", value = Sys.Date(), format = "mm/dd"),
                            timeInput("time1", "Time of First Level (hh:mm):", value = strptime("12:00", format="%H:%M"), seconds = FALSE),
                            numericInput("c1_date", "First Level Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
                            dateInput("date2", "Date of Second Level (t₂):", value = Sys.Date(), format = "mm/dd"),
                            timeInput("time2", "Time of Second Level (hh:mm):", value = strptime("20:00", format="%H:%M"), seconds = FALSE),
                            numericInput("c2_date", "Second Level Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1)
                        ),

                        fluidRow(
                            column(6, numericInput("t_inf", "Infusion Duration (t', hours):", value = 1, min = 0, step = 0.1)),
                            column(6, numericInput("tau", "Dosing Interval (τ, hours):", value = 24, min = 1, step = 0.1))
                        ),
                        numericInput("dose", "Total Dose (mg):", value = 1000, min = 1, step = 1),
                        actionButton("calculate", "Calculate", class = "btn-primary", width = "100%"),
                        actionButton("use_in_opt", "Use results in Optimization", class = "btn-success", width = "100%", style = "margin-top:8px;")
                    )
                ),

                column(8,
                    wellPanel(
                        h4("Calculated PK Parameters", style = "margin-top:0;"),
                        fluidRow(
                            column(6,
                                htmlOutput("ke"),
                                htmlOutput("t_half"),
                                htmlOutput("cmax"),
                                htmlOutput("cmin")
                            ),
                            column(6,
                                htmlOutput("vd"),
                                htmlOutput("clr"),
                                htmlOutput("clr_ml_min"),
                                htmlOutput("auc24")
                            )
                        )
                    ),
                    wellPanel(
                        h4("Concentration-Time Plot", style = "margin-top:0;"),
                        plotOutput("level_plot", height = "360px")
                    )
                )
            )
        ),

        # Dose Optimization tab 
        tabPanel("Dose Optimization",
            fluidRow(
                column(12,
                    p("Use calculated PK parameters (or enter custom values) to optimize dose and predict concentrations.")
                )
            ),
            fluidRow(
                column(4,
                    wellPanel(
                        numericInput("ke_opt", "Elimination Rate (Ke, hr⁻¹):", value = 0.3, min = 0, step = 0.01),
                        numericInput("vd_opt", "Volume of Distribution (Vd, L):", value = 20, min = 0, step = 0.1),
                        numericInput("t_inf_opt", "Infusion Duration (t', hours):", value = 1, min = 0, step = 0.1),
                        numericInput("interval_opt", "Dosing Interval (hours):", value = 24, min = 0, step = 0.1),
                        selectInput("calc_mode", "Calculation Mode:",
                            choices = c("Calculate from AUC Target" = "auc", "Calculate from Dose" = "dose"),
                            selected = "auc"
                        ),
                        conditionalPanel(
                            condition = "input.calc_mode == 'auc'",
                            numericInput("target_auc", "Target AUC24 (mg·h/L):", value = 250, min = 0, step = 1)
                        ),
                        conditionalPanel(
                            condition = "input.calc_mode == 'dose'",
                            numericInput("input_dose", "Dose (mg):", value = 1000, min = 1, step = 1)
                        ),
                        actionButton("optimize_dose", "Calculate", class = "btn-primary", width = "100%")
                    )
                ),
                column(8,
                    wellPanel(
                        h4("Optimization Results", style = "margin-top:0;"),
                        fluidRow(
                            column(6,
                                htmlOutput("opt_dose"),
                                htmlOutput("predicted_peak")
                            ),
                            column(6,
                                htmlOutput("predicted_trough"),
                                htmlOutput("predicted_auc")
                            )
                        )
                    ),
                    wellPanel(
                        plotOutput("opt_plot", height = "320px")
                    )
                )
            )
        )
        ),
    #signature
    br(),
    tags$hr(),
    div("</> Built by Shemual Tsai, last updated 9/15/25. Please reach out to stsai@houstonmethodist.org or any questions or for code.", style = "text-align:center; color:#666; font-size:12px; margin-bottom:10px;")

    )

# Server
server <- function(input, output, session) {
  pk_params <- reactiveValues(ke = NULL, vd = NULL, t_half = NULL, t_inf = NULL, tau = NULL)

  convert_time_to_decimal <- function(time_obj) {
    return(as.numeric(format(time_obj, "%H")) + as.numeric(format(time_obj, "%M")) / 60)
  }

  calculate_time_difference <- function(dose_time, sample_time, infusion_time) {
    dose_time_decimal <- convert_time_to_decimal(dose_time)
    sample_time_decimal <- convert_time_to_decimal(sample_time)
    time_diff <- sample_time_decimal - dose_time_decimal
    return(time_diff)
  }

  # Standard PK Calculator Logic
  observeEvent(input$calculate, {
    dose_datetime <- as.POSIXct(paste(input$dose_date, format(input$dose_time, "%H:%M")))
    datetime1 <- as.POSIXct(paste(input$date1, format(input$time1, "%H:%M")))
    datetime2 <- as.POSIXct(paste(input$date2, format(input$time2, "%H:%M")))

    if (input$input_method == "concentration_time") {
      ke <- log(input$c1 / input$c2) / (input$t2 - input$t1)
      t_half <- 0.693 / ke
      if (ke != 0) {
        cmax <- input$c1 / exp(-ke * (input$t1 - input$t_inf))
      } else {
        cmax <- input$c1
      }

      cmin <- cmax * exp(-ke * (input$tau - input$t_inf))

      vd <- ((input$dose / input$t_inf) * (1 - exp(-ke * input$t_inf))) / (cmax * ke * (1 - exp(-ke * input$tau)))
      clr <- vd * ke
      clr_ml_min <- clr * 16.67
      auc24 <- (input$dose / (ke * vd)) * (24 / input$tau)

      # store PK parameters
      pk_params$ke <- ke
      pk_params$vd <- vd
      pk_params$t_half <- t_half
      pk_params$t_inf <- input$t_inf
      pk_params$tau <- input$tau

      showNotification("PK parameters stored. Click 'Use results in Optimization' to copy them to Dose Optimization tab.", type = "message", duration = 4)

      # Display results 
      output$ke <- renderText({ paste("Elimination Rate Constant (kₑ):", round(ke, 4), "hr⁻¹") })
      output$t_half <- renderText({ paste("Elimination Half-Life (t₁/₂):", round(t_half, 2), "hours") })
      output$cmax <- renderText({ paste("Extrapolated Peak (Cₘₐₓ) at end of infusion (t' = ", input$t_inf, "h):", round(cmax, 2), "mg/L") })
      output$cmin <- renderText({ paste("Extrapolated Trough (Cₘᵢₙ):", round(cmin, 2), "mg/L") })
      output$vd <- renderText({ paste("Volume of Distribution (Vd):", round(vd, 2), "L") })
      output$clr <- renderText({ paste("Clearance (Clᵣ):", round(clr, 2), "L/hr") })
      output$clr_ml_min <- renderText({ paste("Clearance in mL/min:", round(clr_ml_min, 2), "mL/min") })
      output$auc24 <- renderText({ paste("AUC24:", round(auc24, 2), "mg·h/L") })

      # Generate plot
      time_seq <- seq(0, input$tau, by = 0.1)
      conc_seq <- numeric(length(time_seq))
      if (is.null(input$t_inf) || input$t_inf <= 0) {
        time_seq <- seq(0, input$tau, by = 0.1)
        conc_seq <- cmax * exp(-ke * time_seq)
        subtitle_text <- "Infusion time should not be zero or negative"
      } else {
        time_seq <- seq(input$t_inf, input$tau, by = 0.1)
        conc_seq <- cmax * exp(-ke * (time_seq - input$t_inf))
        subtitle_text <- paste0("Curve shown from end of infusion")
      }

      plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
      output$level_plot <- renderPlot({
        ggplot(plot_data, aes(x = Time, y = Concentration)) +
          geom_line(color = "blue", size = 1) +
          geom_point(aes(x = input$t1, y = input$c1), color = "red", size = 3) +
          geom_point(aes(x = input$t2, y = input$c2), color = "red", size = 3) +
          labs(title = "Aminoglycoside PK curve",
               subtitle = subtitle_text,
               x = "Time (hours)", y = "Concentration (mg/L)") +
          theme_minimal()
      })

    } else if (input$input_method == "date_time_dose") {
      t1_diff <- as.numeric(difftime(datetime1, dose_datetime, units = "hours"))
      t2_diff <- as.numeric(difftime(datetime2, dose_datetime, units = "hours"))
      ke <- log(input$c1_date / input$c2_date) / (t2_diff - t1_diff)
      t_half <- 0.693 / ke

      if (ke != 0) {
        cmax <- input$c1_date / exp(-ke * (t1_diff - input$t_inf))
      } else {
        cmax <- input$c1_date
      }

      cmin <- cmax * exp(-ke * (input$tau - input$t_inf))
      vd <- ((input$dose / input$t_inf) * (1 - exp(-ke * input$t_inf))) / (cmax * ke * (1 - exp(-ke * input$tau)))
      clr <- vd * ke
      clr_ml_min <- clr * 16.67
      auc24 <- (input$dose / (ke * vd)) * (24 / input$tau)

      # store PK parameters   
      pk_params$vd <- vd
      pk_params$t_half <- t_half
      pk_params$t_inf <- input$t_inf
      pk_params$tau <- input$tau

      showNotification("PK parameters stored. Click 'Use results in Optimization' to copy them to Dose Optimization tab.", type = "message", duration = 4)

      # Display results
      output$ke <- renderText({ paste("Elimination Rate Constant (kₑ):", round(ke, 4), "hr⁻¹") })
      output$t_half <- renderText({ paste("Elimination Half-Life (t₁/₂):", round(t_half, 2), "hours") })
      output$cmax <- renderText({ paste("Extrapolated Peak (Cₘₐₓ) at end of infusion (t' = ", input$t_inf, "h):", round(cmax, 2), "mg/L") })
      output$cmin <- renderText({ paste("Extrapolated Trough (Cₘᵢₙ):", round(cmin, 2), "mg/L") })
      output$vd <- renderText({ paste("Volume of Distribution (Vd):", round(vd, 2), "L") })
      output$clr <- renderText({ paste("Clearance (Clᵣ):", round(clr, 2), "L/hr") })
      output$clr_ml_min <- renderText({ paste("Clearance in mL/min:", round(clr_ml_min, 2), "mL/min") })
      output$auc24 <- renderText({ paste("AUC24:", round(auc24, 2), "mg·h/L") })

      # Plot 
      time_seq <- seq(0, input$tau, by = 0.1)
      conc_seq <- numeric(length(time_seq))
      if (is.null(input$t_inf) || input$t_inf <= 0) {
        time_seq <- seq(0, input$tau, by = 0.1)
        conc_seq <- cmax * exp(-ke * time_seq)
        subtitle_text <- "Infusion time should not be zero or negative"
      } else {
        time_seq <- seq(input$t_inf, input$tau, by = 0.1)
        conc_seq <- cmax * exp(-ke * (time_seq - input$t_inf))
        subtitle_text <- paste0("Curve shown from end of infusion")
      }

      plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
      output$level_plot <- renderPlot({
        ggplot(plot_data, aes(x = Time, y = Concentration)) +
          geom_line(color = "blue", size = 1) +
          geom_point(aes(x = t1_diff, y = input$c1_date), color = "red", size = 3) +
          geom_point(aes(x = t2_diff, y = input$c2_date), color = "red", size = 3) +
          labs(title = "Aminoglycoside PK Curve",
               subtitle = subtitle_text,
               x = "Time (hours)", y = "Concentration (mg/L)") +
          theme_minimal()
      })
    }
  })

  # Dose Optimization Logic
  observeEvent(input$optimize_dose, {

    calc_peak <- function(dose, vd, ke, t_inf) {
      (dose/(vd * ke * t_inf)) * (1 - exp(-ke * t_inf))
    }

    calc_trough <- function(peak, ke, interval, t_inf) {
      peak * exp(-ke * (interval - t_inf))
    }

    calc_auc24 <- function(dose, vd, ke, interval) {
      (dose / (ke * vd)) * (24 / interval)
    }

    if (input$calc_mode == "peaktrough") {
      new_dose <- (input$vd_opt * input$ke_opt * (1 - exp(-input$ke_opt * input$t_inf_opt)) *
                    input$target_cmax * input$t_inf_opt) / (1 - exp(-input$ke_opt * input$t_inf_opt))

      pred_peak <- input$target_cmax
      pred_trough <- pred_peak * exp(-input$ke_opt * (input$interval_opt - input$t_inf_opt))
      pred_auc24 <- calc_auc24(new_dose, input$vd_opt, input$ke_opt, input$interval_opt)

    } else if (input$calc_mode == "auc") {
      new_dose <- input$target_auc * input$ke_opt * input$vd_opt * (input$interval_opt / 24)

      pred_peak <- calc_peak(new_dose, input$vd_opt, input$ke_opt, input$t_inf_opt)
      pred_trough <- calc_trough(pred_peak, input$ke_opt, input$interval_opt, input$t_inf_opt)
      pred_auc24 <- input$target_auc

    } else { 
      new_dose <- input$input_dose

      pred_peak <- calc_peak(new_dose, input$vd_opt, input$ke_opt, input$t_inf_opt)
      pred_trough <- calc_trough(pred_peak, input$ke_opt, input$interval_opt, input$t_inf_opt)
      pred_auc24 <- calc_auc24(new_dose, input$vd_opt, input$ke_opt, input$interval_opt)
    }

    output$opt_dose <- renderText({ paste("Dose:", round(new_dose, 1), "mg") })
    output$predicted_peak <- renderText({ paste("Predicted Peak:", round(pred_peak, 2), "mg/L") })
    output$predicted_trough <- renderText({ paste("Predicted Trough:", round(pred_trough, 2), "mg/L") })
    output$predicted_auc <- renderText({ paste("Predicted AUC24:", round(pred_auc24, 1), "mg·h/L") })

    # Plot
    time_seq <- seq(0, input$interval_opt, by = 0.1)
    conc_seq <- numeric(length(time_seq))

    if (is.null(input$t_inf_opt) || input$t_inf_opt <= 0) {
      conc_seq <- pred_peak * exp(-input$ke_opt * time_seq)
    } else {
      infusion_idx <- which(time_seq <= input$t_inf_opt)
      post_idx <- which(time_seq > input$t_inf_opt)

      if (length(infusion_idx) > 0) {
      conc_seq[infusion_idx] <- pred_peak * (time_seq[infusion_idx] / input$t_inf_opt)
      }
      if (length(post_idx) > 0) {
      conc_seq[post_idx] <- pred_peak * exp(-input$ke_opt * (time_seq[post_idx] - input$t_inf_opt))
      }
      end_idx <- which(abs(time_seq - input$t_inf_opt) < 1e-8)
      if (length(end_idx) == 0) {
      nearest <- which.min(abs(time_seq - input$t_inf_opt))
      conc_seq[nearest] <- pred_peak
      } else {
      conc_seq[end_idx] <- pred_peak
      }
    }

    plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
    output$opt_plot <- renderPlot({
      ggplot(plot_data, aes(x = Time, y = Concentration)) +
      geom_line(color = "blue", size = 1) +
      geom_hline(yintercept = pred_peak, linetype = "dashed", color = "red", alpha = 0.5) +
      geom_hline(yintercept = pred_trough, linetype = "dashed", color = "red", alpha = 0.5) +
      labs(title = "Predicted Concentration-Time Curve",
         x = "Time (hours)", y = "Concentration (mg/L)") +
      theme_minimal()
    })
    })

  # Copy calculated PK parameters into Dose Optimization tab
  observeEvent(input$use_in_opt, {
    req(pk_params$ke, pk_params$vd)  
    updateNumericInput(session, "ke_opt", value = round(pk_params$ke, 4))
    updateNumericInput(session, "vd_opt", value = round(pk_params$vd, 2))
    if (!is.null(pk_params$t_inf)) updateNumericInput(session, "t_inf_opt", value = pk_params$t_inf)
    if (!is.null(pk_params$tau)) updateNumericInput(session, "interval_opt", value = pk_params$tau)
    showNotification("Patient-specific parameters copied to Dose Optimization tab", type = "message", duration = 3)
  })
}


# Run the application
shinyApp(ui, server)
