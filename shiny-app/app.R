library(shiny)

library(bslib)

library(ggplot2)

library(tidyr)

source("server-plots.R")

ui <- page_sidebar(
  
  title = "Decoherence in a Hadamard Quantum Random Walk on a Line",
  
  sidebar = sidebar(
    
    width = "360px",
    
    open = "desktop",
    
    accordion(
      
      multiple = TRUE,
      
      accordion_panel(
        
        "Common settings",
        
        numericInput("T", "Number of steps (T):", value = 20, min = 1, step = 1),
        
        selectInput(
          
          "init_coin",
          
          "Initial coin state:",
          
          choices = c("symmetric", "right", "left"),
          
          selected = "symmetric"
        )
      ),
      
      accordion_panel(
        
        "Noisy QRW settings",
        
        checkboxInput("show_noisy", "Display noisy QRW", value = TRUE),
        
        numericInput("N", "Number of noisy trajectories (N):", value = 200, min = 1, step = 1),
        
        selectInput(
          
          "channel",
          
          "Decoherence channel:",
          
          choices = c("dephasing", "depolarizing"),
          
          selected = "dephasing"
        ),
        
        sliderInput(
          
          "p",
          
          "Decoherence strength (p):",
          
          min = 0, max = 1, value = 0.1, step = 0.01
        ),
        
        numericInput("seed", "Seed:", value = 1, min = 1, step = 1)
      ),
      
      accordion_panel(
        
        "Noiseless QRW settings",
        
        checkboxInput("show_noiseless", "Display noiseless QRW", value = TRUE),
        
        helpText("Uses the common initial coin state selected above.")
      ),
      
      accordion_panel(
        
        "Classical RW settings",
        
        checkboxInput("show_classical", "Display classical RW", value = TRUE),
        
        helpText("The classical symmetric random walk uses the same number of steps T.")
      ),
      
      accordion_panel(
        
        "Display settings",
        
        radioButtons(
          
          "display_mode",
          
          "Display mode:",
          
          choices = c("Overlay" = "overlay", "Separate graphs" = "separate"),
          
          selected = "overlay"
        ),
        
        helpText("The data table always matches the currently active graph tab.")
      )
    ),
    
    actionButton("run", "Run simulation", class = "btn-primary")
  ),
  
  navset_card_tab(
    
    id = "main_tabs",
    
    nav_panel(
      
      "Final Distribution",
      
      conditionalPanel(
        
        condition = "input.display_mode == 'overlay'",
        
        plotOutput("dist_overlay_plot", height = "420px")
      ),
      
      conditionalPanel(
        
        condition = "input.display_mode == 'separate'",
        
        uiOutput("dist_separate_ui")
      )
    ),
    
    nav_panel(
      
      "Variance vs Time",
      
      conditionalPanel(
        
        condition = "input.display_mode == 'overlay'",
        
        plotOutput("var_overlay_plot", height = "420px")
      ),
      
      conditionalPanel(
        
        condition = "input.display_mode == 'separate'",
        
        uiOutput("var_separate_ui")
      )
    ),
    
    nav_panel(
      
      "Displayed Data",
      
      h4(textOutput("active_table_title")),
      
      tableOutput("active_table")
    ),
    
    nav_panel(
      
      "Final Summary",
      
      h4("Final Summary Table"),
      
      tableOutput("summary_tbl")
    )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    noisy <- sim_noisy_qrw(
      
      T = input$T,
      
      N = input$N,
      
      channel = input$channel,
      
      init_coin = input$init_coin,
      
      p = input$p,
      
      seed = input$seed
    )
    
    noiseless <- sim_noiseless_qrw(
      
      T = input$T,
      
      init_coin = input$init_coin
    )
    
    classical <- sim_srw(input$T)
    
    list(
      
      noisy = noisy,
      
      noiseless = noiseless,
      
      classical = classical,
      
      dists = build_dist(noisy, noiseless, classical),
      
      vars = build_variance_overlay(noisy, noiseless, classical),
      
      summary_tbl = build_summary_table(noisy, noiseless, classical)
    )
    
  }, ignoreNULL = FALSE)
  
  selected_models <- reactive({
    
    mods <- character(0)
    
    if (isTRUE(input$show_noisy)) mods <- c(mods, "noisy")
    
    if (isTRUE(input$show_noiseless)) mods <- c(mods, "noiseless")
    
    if (isTRUE(input$show_classical)) mods <- c(mods, "classical")
    
    req(length(mods) > 0)
    
    mods
  })
  
  dist_plot_data <- reactive({
    
    res <- results()
    
    req(res)
    
    make_dist_plot_data(res$dists, selected_models())
  })
  
  var_plot_data <- reactive({
    
    res <- results()
    
    req(res)
    
    make_var_plot_data(res$vars, selected_models())
  })
  
  output$dist_overlay_plot <- renderPlot({
    
    plot_dist_overlay(
      
      plot_df = dist_plot_data(),
      
      selected_models = selected_models(),
      
      T = input$T,
      
      N = input$N,
      
      channel = input$channel,
      
      p = input$p
    )
  })
  
  output$var_overlay_plot <- renderPlot({
    
    plot_var_overlay(
      
      plot_df = var_plot_data(),
      
      selected_models = selected_models()
    )
  })
  
  output$dist_separate_ui <- renderUI({
    
    req(length(selected_models()) > 0)
    
    plot_list <- lapply(selected_models(), function(m) {
      
      plotOutput(paste0("dist_", m), height = "350px")
    })
    
    do.call(tagList, plot_list)
  })
  
  output$var_separate_ui <- renderUI({
    
    req(length(selected_models()) > 0)
    
    plot_list <- lapply(selected_models(), function(m) {
      
      plotOutput(paste0("var_", m), height = "350px")
    })
    
    do.call(tagList, plot_list)
  })
  
  observe({
    
    res <- results()
    
    req(res)
    
    for (m in c("noisy", "noiseless", "classical")) {
      
      local({
        
        model_name <- m
        
        output[[paste0("dist_", model_name)]] <- renderPlot({
          
          req(model_name %in% selected_models())
          
          plot_dist_single(res$dists, model_name, input$T, input$init_coin)
        })
        
        output[[paste0("var_", model_name)]] <- renderPlot({
          
          req(model_name %in% selected_models())
          
          plot_var_single(res$vars, model_name)
        })
      })
    }
  })
  
  output$active_table_title <- renderText({
    
    if (input$main_tabs == "Final Distribution") {
      
      "Displayed Data Table: Final Distribution"
      
    } else if (input$main_tabs == "Variance vs Time") {
      
      "Displayed Data Table: Variance vs Time"
      
    } else {
      
      "Displayed Data Table"
    }
  })
  
  output$active_table <- renderTable({
    
    res <- results()
    
    req(res)
    
    req(length(selected_models()) > 0)
    
    make_active_table(
      
      results = res,
      
      main_tab = input$main_tabs,
      
      selected_models = selected_models()
    )
    
  }, rownames = FALSE)
  
  output$summary_tbl <- renderTable({
    
    res <- results()
    
    req(res)
    
    res$summary_tbl
    
  }, rownames = FALSE)
}

shinyApp(ui, server)