# Set CRAN mirror
options(repos = c(CRAN = "http://cran.rstudio.com"))

# Install required packages if not installed
required_packages <- c("shiny", "ggplot2", "shinydashboard", "shinyTime")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
library(shiny)
library(ggplot2)
library(shinydashboard)
library(shinyTime) 

# Set up UI
ui <- dashboardPage(
  dashboardHeader(title = "Aminoglycoside Calculator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Standard PK Calculator", tabName = "standard", icon = icon("calculator")),
      menuItem("Dose Optimization", tabName = "optimize", icon = icon("sliders"))
    )
  ),
  dashboardBody(
    tabItems(
      # Standard PK Calculator tab
      tabItem(tabName = "standard",
              fluidRow(
                box(
                  title = "Input Parameters",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 4,
                  radioButtons("input_method", "Select Input Method:",
                               choices = c("Enter Concentrations and Times" = "concentration_time",
                                           "Enter Date/Time and Dose" = "date_time_dose")),
                  conditionalPanel(
                    condition = "input.input_method == 'concentration_time'",
                    numericInput("c1", "First Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
                    numericInput("c2", "Second Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1),
                    numericInput("t1", "Time of First Level (t₁, hours) after END of infusion:", value = 1, min = 0, step = 0.1),
                    numericInput("t2", "Time of Second Level (t₂, hours) after END of infusion:", value = 8, min = 0, step = 0.1)
                  ),
                  conditionalPanel(
                    condition = "input.input_method == 'date_time_dose'",
                    dateInput("dose_date", "Date of Dose Administration:", value = Sys.Date()),
                    timeInput("dose_time", "Time of Dose Administration:", value = strptime("12:00", format="%H:%M")),
                    dateInput("date1", "Date of First Level (t₁):", value = Sys.Date()),
                    timeInput("time1", "Time of First Level (t₁):", value = strptime("13:00", format="%H:%M")),
                    numericInput("c1_date", "First Level Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
                    dateInput("date2", "Date of Second Level (t₂):", value = Sys.Date()),
                    timeInput("time2", "Time of Second Level (t₂):", value = strptime("20:00", format="%H:%M")),
                    numericInput("c2_date", "Second Level Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1),
                  ),
                  numericInput("t_inf", "Infusion Duration (t', hours):", value = 1, min = 0, step = 0.1),
                  numericInput("tau", "Dosing Interval (τ, hours):", value = 24, min = 1, step = 0.1),
                    numericInput("dose", "Total Dose (mg):", value = 400, min = 1, step = 1),
                  actionButton("calculate", "Calculate", class = "btn-primary")
                ),
                box(
                  title = "Results",
                  status = "info",
                  solidHeader = TRUE,
                  width = 8,
                  verbatimTextOutput("ke"),
                  verbatimTextOutput("t_half"),
                  verbatimTextOutput("cmax"),
                  verbatimTextOutput("cmin"),
                  verbatimTextOutput("vd"),
                  verbatimTextOutput("clr"),
                  verbatimTextOutput("clr_ml_min"),
                  verbatimTextOutput("auc24")
                )
              ),
              fluidRow(
                box(
                  title = "Concentration-Time Curve",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  plotOutput("level_plot")
                )
              )
      ),
      
      # Dose Optimization tab
      tabItem(tabName = "optimize",
              fluidRow(
                box(
                  title = "PK Parameters",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 4,
                  numericInput("ke_opt", "Elimination Rate (Ke, hr⁻¹):", value = 0.2, min = 0, step = 0.01),
                  numericInput("vd_opt", "Volume of Distribution (Vd, L):", value = 20, min = 0, step = 0.1),
                  numericInput("t_inf_opt", "Infusion Duration (t', hours):", value = 1, min = 0, step = 0.1),
                  numericInput("interval_opt", "Dosing Interval (hours):", value = 24, min = 0, step = 0.1),
                  hr(),
                  selectInput("calc_mode", "Calculation Mode:",
                              choices = c("Calculate from AUC Target" = "auc",
                                          "Calculate from Dose" = "dose")),
                  conditionalPanel(
                    condition = "input.calc_mode == 'auc'",
                    numericInput("target_auc", "Target AUC24 (mg·h/L):", value = 400, min = 0, step = 1)
                  ),
                  conditionalPanel(
                    condition = "input.calc_mode == 'dose'",
                    numericInput("input_dose", "Dose (mg):", value = 400, min = 1, step = 1)
                  ),
                  actionButton("optimize_dose", "Calculate", class = "btn-primary")
                ),
                box(
                  title = "Results",
                  status = "info",
                  solidHeader = TRUE,
                  width = 8,
                  verbatimTextOutput("opt_dose"),
                  verbatimTextOutput("predicted_peak"),
                  verbatimTextOutput("predicted_trough"),
                  verbatimTextOutput("predicted_auc")
                )
              ),
              fluidRow(
                box(
                  title = "Predicted Concentration-Time Curve",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  plotOutput("opt_plot")
                )
              )
      )
    )
  )
)

# Set up server
server <- function(input, output, session) {
  # Function to convert time from POSIXct (HH:MM) to decimal time (fraction of a day)
  convert_time_to_decimal <- function(time_obj) {
    return(as.numeric(format(time_obj, "%H")) + as.numeric(format(time_obj, "%M")) / 60)
  }
  # Function to calculate the time differences
  calculate_time_difference <- function(dose_time, sample_time, infusion_time) {
    # Convert to numeric time (decimal hours)
    dose_time_decimal <- convert_time_to_decimal(dose_time)
    sample_time_decimal <- convert_time_to_decimal(sample_time)
    time_diff <- sample_time_decimal - dose_time_decimal
    return(time_diff)
  }
  
  # PK Calculator Logic
  observeEvent(input$calculate, {
    # Convert date and time inputs to POSIXct objects
    dose_datetime <- as.POSIXct(paste(input$dose_date, format(input$dose_time, "%H:%M:%S")))
    datetime1 <- as.POSIXct(paste(input$date1, format(input$time1, "%H:%M:%S")))
    datetime2 <- as.POSIXct(paste(input$date2, format(input$time2, "%H:%M:%S")))
    
    if (input$input_method == "concentration_time") {
      # Original PK calculations based on concentrations and times
      ke <- log(input$c1 / input$c2) / (input$t2 - input$t1)
      t_half <- 0.693 / ke
  
      if (ke != 0) {
        cmax <- input$c1 / exp(-ke * (input$t1 - input$t_inf))
      } else {
        cmax <- input$c1  # if ke is zero, cmax equals c1
      }
      
      cmin <- cmax * exp(-ke * (input$tau - input$t_inf))
      vd <- input$dose / cmax
      clr <- vd * ke
      clr_ml_min <- clr * 16.67
      auc24 <- (input$dose / (ke * vd)) * (24 / input$tau)
      
      # Display results
      output$ke <- renderText({ paste("Elimination Rate Constant (kₑ):", round(ke, 4), "hr⁻¹") })
      output$t_half <- renderText({ paste("Elimination Half-Life (t₁/₂):", round(t_half, 2), "hours") })
      output$cmax <- renderText({ paste("Extrapolated Peak (Cₘₐₓ):", round(cmax, 2), "mg/L") })
      output$cmin <- renderText({ paste("Extrapolated Trough (Cₘᵢₙ):", round(cmin, 2), "mg/L") })
      output$vd <- renderText({ paste("Volume of Distribution (Vd):", round(vd, 2), "L") })
      output$clr <- renderText({ paste("Clearance (Clᵣ):", round(clr, 2), "L/hr") })
      output$clr_ml_min <- renderText({ paste("Clearance in mL/min:", round(clr_ml_min, 2), "mL/min") })
      output$auc24 <- renderText({ paste("AUC24:", round(auc24, 2), "mg·h/L") })
      
      # Generate plot
      time_seq <- seq(0, input$tau, by = 0.1)
      conc_seq <- cmax * exp(-ke * (time_seq))
      
      plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
      output$level_plot <- renderPlot({
        ggplot(plot_data, aes(x = Time, y = Concentration)) +
          geom_line(color = "blue", size = 1) +
       geom_point(aes(x = input$t1-input$t_inf, y = input$c1), color = "red", size = 3) +
        geom_point(aes(x = input$t2-input$t_inf, y = input$c2), color = "red", size = 3) +
          labs(title = "Aminoglycoside Concentration-Time Curve",
               x = "Time (hours)", y = "Concentration (mg/L)") +
          theme_minimal()
      })
    } else if (input$input_method == "date_time_dose") {
      # Calculate time differences in hours
      t1_diff <- as.numeric(difftime(datetime1, dose_datetime, units = "hours")) - input$t_inf
      t2_diff <- as.numeric(difftime(datetime2, dose_datetime, units = "hours")) - input$t_inf
    
      ke <- log(input$c1_date / input$c2_date) / (t2_diff - t1_diff)
      t_half <- 0.693 / ke
      
      if (ke != 0) {
        cmax <- input$c1_date / exp(-ke * (t1_diff))
      } else {
        cmax <- input$c1_date  # if ke is zero, cmax equals c1
      }
      
      cmin <- cmax * exp(-ke * (input$tau - input$t_inf))
      vd <- input$dose / cmax
      clr <- vd * ke
      clr_ml_min <- clr * 16.67
      auc24 <- (input$dose / (ke * vd)) * (24 / input$tau)
      
      # Display results
      output$ke <- renderText({ paste("Elimination Rate Constant (kₑ):", round(ke, 4), "hr⁻¹") })
      output$t_half <- renderText({ paste("Elimination Half-Life (t₁/₂):", round(t_half, 2), "hours") })
      output$cmax <- renderText({ paste("Extrapolated Peak (Cₘₐₓ):", round(cmax, 2), "mg/L") })
      output$cmin <- renderText({ paste("Extrapolated Trough (Cₘᵢₙ):", round(cmin, 2), "mg/L") })
      output$vd <- renderText({ paste("Volume of Distribution (Vd):", round(vd, 2), "L") })
      output$clr <- renderText({ paste("Clearance (Clᵣ):", round(clr, 2), "L/hr") })
      output$clr_ml_min <- renderText({ paste("Clearance in mL/min:", round(clr_ml_min, 2), "mL/min") })
      output$auc24 <- renderText({ paste("AUC24:", round(auc24, 2), "mg·h/L") })
      
      # Generate plot
      time_seq <- seq(0, input$tau, by = 0.1)
      conc_seq <- cmax * exp(-ke * (time_seq))
      
      plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
      output$level_plot <- renderPlot({
        ggplot(plot_data, aes(x = Time, y = Concentration)) +
          geom_line(color = "blue", size = 1) +
          geom_point(aes(x = t1_diff, y = input$c1_date), color = "red", size = 3) +
          geom_point(aes(x = t2_diff, y = input$c2_date), color = "red", size = 3) +
          labs(title = "Aminoglycoside Concentration-Time Curve",
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
    
    # Calculate based on selected mode
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
    
    # Display optimization results
    output$opt_dose <- renderText({ paste("Dose:", round(new_dose, 1), "mg") })
    output$predicted_peak <- renderText({ paste("Predicted Peak:", round(pred_peak, 2), "mg/L") })
    output$predicted_trough <- renderText({ paste("Predicted Trough:", round(pred_trough, 2), "mg/L") })
    output$predicted_auc <- renderText({ paste("Predicted AUC24:", round(pred_auc24, 1), "mg·h/L") })
    
    # Generate optimization plot
    time_seq <- seq(0, input$interval_opt, by = 0.1)
    conc_seq <- pred_peak * exp(-input$ke_opt * (time_seq - input$t_inf_opt))
    conc_seq[time_seq < input$t_inf_opt] <- 0  # Zero concentration during infusion
    
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
}

shinyApp(ui, server)
