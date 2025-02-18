# Set CRAN mirror
options(repos = c(CRAN = "http://cran.rstudio.com"))

# Install required packages if not installed
required_packages <- c("shiny", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
library(shiny)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("Aminoglycoside AUC & Pharmacokinetics Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("c1", "First Aminoglycoside Concentration (C₁, mg/L):", value = 8, min = 0, step = 0.1),
      numericInput("c2", "Second Aminoglycoside Concentration (C₂, mg/L):", value = 1, min = 0, step = 0.1),
      numericInput("t1", "Time of First Level (t₁, hours):", value = 1, min = 0, step = 0.1),
      numericInput("t2", "Time of Second Level (t₂, hours):", value = 8, min = 0, step = 0.1),
      numericInput("t_inf", "Infusion Duration (t’, hours):", value = 0.5, min = 0, step = 0.1),
      numericInput("tau", "Dosing Interval (τ, hours):", value = 12, min = 1, step = 0.1),
      numericInput("dose", "Total Dose (mg):", value = 400, min = 1, step = 1),
      actionButton("calculate", "Calculate")
    ),
    
    mainPanel(
      h3("Results"),
      verbatimTextOutput("ke"),
      verbatimTextOutput("t_half"),
      verbatimTextOutput("cmax"),
      verbatimTextOutput("cmin"),
      verbatimTextOutput("vd"),
      verbatimTextOutput("clr"),
      verbatimTextOutput("clr_ml_min"),
      verbatimTextOutput("auc24"),
      plotOutput("level_plot")
    )
  )
)

# Define server
server <- function(input, output, session) {
  observeEvent(input$calculate, {
    # Perform pharmacokinetic calculations
    ke <- log(input$c1 / input$c2) / (input$t2 - input$t1)
    t_half <- 0.693 / ke
    cmax <- input$c1 / exp(-ke * (input$t1 - input$t_inf))
    cmin <- cmax * exp(-ke * (input$tau - input$t_inf))
    vd <- input$dose / cmax
    clr <- vd * ke
    clr_ml_min <- clr * 16.67
    auc24 <- (input$dose / (ke * vd)) * (24 / input$tau)
    
    # Display the results
    output$ke <- renderText({ paste("Elimination Rate Constant (kₑ):", round(ke, 4), "1/hr") })
    output$t_half <- renderText({ paste("Elimination Half-Life (t₁/₂):", round(t_half, 2), "hours") })
    output$cmax <- renderText({ paste("Extrapolated Peak (Cₘₐₓ):", round(cmax, 2), "mg/L") })
    output$cmin <- renderText({ paste("Extrapolated Trough (Cₘᵢₙ):", round(cmin, 2), "mg/L") })
    output$vd <- renderText({ paste("Volume of Distribution (Vd):", round(vd, 2), "L") })
    output$clr <- renderText({ paste("Clearance (Clᵣ):", round(clr, 2), "L/hr") })
    output$clr_ml_min <- renderText({ paste("Clearance in mL/min:", round(clr_ml_min, 2), "mL/min") })
    output$auc24 <- renderText({ paste("AUC24:", round(auc24, 2), "mg·h/L") })
    
    # Data for plotting
    time_seq <- seq(0, input$tau, by = 0.1)
    conc_seq <- cmax * exp(-ke * (time_seq))
    
    # Generate concentration-time curve
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
}

# Run the app
shinyApp(ui = ui, server = server)
