# Set CRAN mirror
options(repos = c(CRAN = "http://cran.rstudio.com"))

# Install required packages if not installed
required_packages <- c("shiny", "ggplot2", "shinydashboard")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
library(shiny)
library(ggplot2)
library(shinydashboard)

# Define UI
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
            numericInput("c1", "First Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
            numericInput("c2", "Second Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1),
            numericInput("t1", "Time of First Level (t₁, hours):", value = 1, min = 0, step = 0.1),
            numericInput("t2", "Time of Second Level (t₂, hours):", value = 8, min = 0, step = 0.1),
            numericInput("t_inf", "Infusion Duration (t', hours):", value = 0.5, min = 0, step = 0.1),
            numericInput("tau", "Dosing Interval (τ, hours):", value = 12, min = 1, step = 0.1),
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
            numericInput("t_inf_opt", "Infusion Duration (t', hours):", value = 0.5, min = 0, step = 0.1),
            numericInput("interval_opt", "Dosing Interval (hours):", value = 24, min = 0, step = 0.1),
            hr(),
            selectInput("calc_mode", "Calculation Mode:",
                       choices = c("Calculate from Peak/Trough Targets" = "peaktrough",
                                 "Calculate from AUC Target" = "auc",
                                 "Calculate from Dose" = "dose")),
            conditionalPanel(
              condition = "input.calc_mode == 'peaktrough'",
              numericInput("target_cmax", "Target Peak (Cmax, mg/L):", value = 20, min = 0, step = 0.1),
              numericInput("target_cmin", "Target Trough (Cmin, mg/L):", value = 1, min = 0, step = 0.1)
            ),
            conditionalPanel(
              condition = "input.calc_mode == 'auc'",
              numericInput("target_auc", "Target AUC24 (mg·h/L):", value = 400, min = 0, step = 1)
            ),
            conditionalPanel(
              condition = "input.calc_mode == 'dose'",
              numericInput("input_dose", "Dose (mg):", value = 400, min = 0, step = 1)
            ),
            actionButton("optimize_dose", "Calculate", class = "btn-primary")
          ),
          box(
            title = "Optimization Results",
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

# Define server
server <- function(input, output, session) {
  # Standard PK Calculator Logic
  observeEvent(input$calculate, {
    # Original PK calculations
    ke <- log(input$c1 / input$c2) / (input$t2 - input$t1)
    t_half <- 0.693 / ke
    cmax <- input$c1 / exp(-ke * (input$t1 - input$t_inf))
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
        geom_point(aes(x = input$t1, y = input$c1), color = "red", size = 3) +
        geom_point(aes(x = input$t2, y = input$c2), color = "red", size = 3) +
        labs(title = "Aminoglycoside Concentration-Time Curve",
             x = "Time (hours)", y = "Concentration (mg/L)") +
        theme_minimal()
    })
  })
  
  # Dose Optimization Logic
  observeEvent(input$optimize_dose, {
    # Function to calculate peak from dose
    calc_peak <- function(dose, vd, ke, t_inf) {
      (dose/(vd * ke * t_inf)) * (1 - exp(-ke * t_inf))
    }
    
    # Function to calculate trough from peak
    calc_trough <- function(peak, ke, interval, t_inf) {
      peak * exp(-ke * (interval - t_inf))
    }
    
    # Function to calculate AUC24 from dose
    calc_auc24 <- function(dose, vd, ke, interval) {
      (dose / (ke * vd)) * (24 / interval)
    }
    
    if (input$calc_mode == "peaktrough") {
      # Calculate dose from target peak
      new_dose <- (input$vd_opt * input$ke_opt * (1 - exp(-input$ke_opt * input$t_inf_opt)) * 
                   input$target_cmax * input$t_inf_opt) / (1 - exp(-input$ke_opt * input$t_inf_opt))
      
      pred_peak <- input$target_cmax
      pred_trough <- pred_peak * exp(-input$ke_opt * (input$interval_opt - input$t_inf_opt))
      pred_auc24 <- calc_auc24(new_dose, input$vd_opt, input$ke_opt, input$interval_opt)
      
    } else if (input$calc_mode == "auc") {
      # Calculate dose from target AUC24
      new_dose <- input$target_auc * input$ke_opt * input$vd_opt * (input$interval_opt / 24)
      
      pred_peak <- calc_peak(new_dose, input$vd_opt, input$ke_opt, input$t_inf_opt)
      pred_trough <- calc_trough(pred_peak, input$ke_opt, input$interval_opt, input$t_inf_opt)
      pred_auc24 <- input$target_auc
      
    } else { # dose mode
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
    time_seq <- seq(0, new_interval, by = 0.1)
    conc_seq <- pred_peak * exp(-input$ke_opt * (time_seq - input$t_inf_opt))
    
    plot_data <- data.frame(Time = time_seq, Concentration = conc_seq)
    output$opt_plot <- renderPlot({
      ggplot(plot_data, aes(x = Time, y = Concentration)) +
        geom_line(color = "blue", size = 1) +
        geom_hline(yintercept = input$target_cmax, linetype = "dashed", color = "red") +
        geom_hline(yintercept = input$target_cmin, linetype = "dashed", color = "red") +
        labs(title = "Predicted Concentration-Time Curve with Optimal Dosing",
             x = "Time (hours)", y = "Concentration (mg/L)") +
        theme_minimal()
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
