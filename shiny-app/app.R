library(shiny)

library(bslib)

library(ggplot2)

library(tidyr)

library(plotly)

source("server-plots.R")

ui <- page_sidebar(
  
  title = div(
    
    class = "d-flex justify-content-between align-items-center w-100",
    
    "Decoherence in a Hadamard Quantum Random Walk on a Line",
    
    input_dark_mode(id = "dark")
  ),
  
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
        
        plotlyOutput("dist_overlay_plot", height = "420px")
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
        
        plotlyOutput("var_overlay_plot", height = "420px")
      ),
      
      conditionalPanel(
        
        condition = "input.display_mode == 'separate'",
        
        uiOutput("var_separate_ui")
      )
    ),
    
    nav_panel(
      
      "Final Summary",
      
      h4("Final Summary Table"),
      
      tableOutput("summary_tbl")
    ),
    
    nav_panel(
      "Model Guide",
      
      layout_column_wrap(
        width = 1,
        
        card(
          card_header("Overview"),
          card_body(
            markdown("
            
            This app compares three random-walk models on a one-dimensional line:
            
            - **Noiseless Hadamard quantum random walk (QRW)**  
            - **Noisy Hadamard quantum random walk**, where decoherence is added through a Kraus channel  
            - **Classical symmetric random walk (SRW)**  
            
            The purpose is to show how decoherence changes the behaviour of the quantum walk and can drive it toward more classical behaviour.
            ")
          )
        ),
        
        accordion(
          
          accordion_panel(
            "1. Noiseless quantum random walk",
            withMathJax(
            markdown("
            The noiseless quantum random walk is the fully coherent version of the model.  
            At each step, the walker first undergoes a **Hadamard coin flip** and then a **conditional shift**.
            
            The update is
            
            $$
            U = S(H \\otimes I_p),
            $$
            
            where:
            
            - \\\\(H\\\\)  is the Hadamard coin operator,
            - \\\\(S\\\\)  is the conditional shift,
            - \\\\(I_p\\\\) is the identity on position space.
            
            Because the state evolves coherently, interference occurs between different paths.  
            This usually creates a distribution with strong oscillations and faster spreading than a classical walk.
            ")
            )
          ),
          
          accordion_panel(
            "2. Noisy quantum random walk",
            markdown("
            The noisy quantum random walk adds **decoherence** to the coin state.
            
            Instead of evolving only with the unitary step, the walk also applies a noise channel described by **Kraus operators**.  
            This app includes:
            
            - **Dephasing**: reduces phase coherence while preserving basis populations
            - **Depolarizing**: mixes the coin state more strongly and tends to destroy coherence faster
            
            The noisy model is simulated using **Monte Carlo quantum trajectories** and then averaged across trajectories.
            
            As the decoherence strength **p** increases, the walk typically becomes less quantum and more classical in appearance.
            ")
          ),
          
          accordion_panel(
            "3. Classical symmetric random walk",
            withMathJax(
            markdown("
            
            The classical symmetric random walk is the baseline comparison model.
            
            At each step, the walker moves:
            
            - one unit left with probability $$p=\\frac{1}{2}$$
            - one unit right with probability $$p=\\frac{1}{2}$$
            
            Unlike the quantum walk, there is no superposition and no interference.  
            Its distribution is smoother and spreads more slowly than the noiseless quantum walk.
            ")
            )
          ),
                      
          accordion_panel(
            "4. How to read the final distribution graph",
            markdown("
            The **final distribution graph** shows the probability of the walker being at each position after the chosen number of steps.
            
            How to interpret it:
            
            - **Wide spread** means the walk has explored positions farther from the origin
            - **Sharp peaks and oscillations** are typical signatures of quantum interference
            - **Smoother bell shaped structure** is more classical in character
            - If the noisy walk looks more like the classical walk as **p** increases, that suggests a quantum to classical transition
            
            Typical behaviour:
            
            - **Noiseless QRW**: oscillatory and widely spread
            - **Noisy QRW**: intermediate behaviour depending on **p**
            - **Classical RW**: smoother and more concentrated near the centre
            ")
          ),
          
          accordion_panel(
            "5. How to read the variance graph",
            withMathJax(
            markdown("
            
            The **variance vs. time** graph shows how quickly the walk spreads over time.
            
            This is one of the most important summaries of behaviour.
            
            - **Classical random walk** typically has diffusive scaling:
              
            $$
            \\mathrm{Var}(X_t) \\sim t
            $$
            
            - **Noiseless quantum walk** typically has ballistic scaling:
              
            $$
            \\mathrm{Var}(X_t) \\sim t^2
            $$
            
            So:
            
            - a curve growing roughly linearly is more classical
            - a curve growing much faster is more quantum
            - if increasing decoherence causes the noisy curve to move closer to the classical one, that indicates loss of coherence
            ")
            )
          ),
          
          accordion_panel(
            
            "6. Meaning of the controls",
            withMathJax(
            markdown("
            
            ### `T`
            Number of time steps in the walk.
            
            ### `N`
            Number of Monte Carlo trajectories used for the noisy quantum walk.  
            Larger `N` usually gives a more stable average but takes longer to compute.
            
            ### `p`
            Decoherence strength.
            
            - `p = 0` means no decoherence
            - larger `p` means stronger noise
            - stronger noise generally suppresses quantum interference
            
            ### Channel
            The selected noise mechanism for the noisy quantum walk.
            
            - **Dephasing** mainly destroys phase information
            - **Depolarizing** introduces stronger mixing
            
            ### Initial coin state
            The initial state of the quantum coin.  
            This can affect symmetry and the detailed shape of the distribution.
            ")
            )
          ),
          
          accordion_panel(
            "7. General interpretation tips",
            markdown("
            
            When comparing outputs, ask:
            
            1. **How wide is the distribution?**  
               Wider spreading usually indicates stronger quantum behaviour.
            
            2. **Are there oscillations or interference patterns?**  
               These are signatures of coherence.
            
            3. **Does the noisy walk move closer to the classical walk as **p** increases?**  
               If yes, this supports the quantum to classical transition.
            
            4. **How does the variance curve grow over time?**  
               Faster growth suggests ballistic quantum spreading, while slower linear growth suggests classical diffusion.
            ")
          )
        )
      )
    )
    
    
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    withProgress(message = "Running simulation...", value = 0, expr = {
      
      
      incProgress(0.1, detail = "Simulating noiseless QRW...")
      
      noiseless <- sim_noiseless_qrw(
        
        T = input$T,
        
        init_coin = input$init_coin
      )
      
      incProgress(0.4, detail = "Simulating noisy QRW...")
      
      noisy <- sim_noisy_qrw(
        
        T = input$T,
        
        N = input$N,
        
        channel = input$channel,
        
        init_coin = input$init_coin,
        
        p = input$p,
        
        seed = input$seed
      )
      
      incProgress(0.3, detail = "Simulating classical RW...")
      
      classical <- sim_srw(input$T)
      
      incProgress(0.2, detail = "Building results...")
      
      list(
        
        noisy = noisy,
        
        noiseless = noiseless,
        
        classical = classical,
        
        dists = build_dist(noisy, noiseless, classical),
        
        vars = build_variance_overlay(noisy, noiseless, classical),
        
        summary_tbl = build_summary_table(noisy, noiseless, classical)
      )
    })
    
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
  
  output$dist_overlay_plot <- renderPlotly({
    
    dark <- isTRUE(input$dark == "dark")
    
    plot_dist_overlay(
      
      plot_df = dist_plot_data(),
      
      selected_models = selected_models(),
      
      T = input$T,
      
      N = input$N,
      
      channel = input$channel,
      
      p = input$p,
      
      dark = dark
    )
    
  })
  
  output$var_overlay_plot <- renderPlotly({
    
    dark <- isTRUE(input$dark == "dark")
    
    plot_var_overlay(
      
      plot_df = var_plot_data(),
      
      selected_models = selected_models(),
      
      dark = dark
    )
    
  })
  
  output$dist_separate_ui <- renderUI({
    
    req(length(selected_models()) > 0)
    
    plot_list <- lapply(selected_models(), function(m) {
      
      plotlyOutput(paste0("dist_", m), height = "350px")
    })
    
    do.call(tagList, plot_list)
  })
  
  output$var_separate_ui <- renderUI({
    
    req(length(selected_models()) > 0)
    
    plot_list <- lapply(selected_models(), function(m) {
      
      plotlyOutput(paste0("var_", m), height = "350px")
    })
    
    do.call(tagList, plot_list)
  })
  
  observe({
    
    res <- results()
    
    req(res)
    
    for (m in c("noisy", "noiseless", "classical")) {
      
      local({
        
        model_name <- m
        
        output[[paste0("dist_", model_name)]] <- renderPlotly({
          
          req(model_name %in% selected_models())
          
          dark <- isTRUE(input$dark == "dark")
          
          plot_dist_single(
            
            res$dists,
            
            model_name,
            
            input$T,
            
            input$init_coin,
            
            dark = dark
          )
          
        })
        
        output[[paste0("var_", model_name)]] <- renderPlotly({
          
          req(model_name %in% selected_models())
          
          dark <- isTRUE(input$dark == "dark")
          
          plot_var_single(
            
            res$vars,
            
            model_name,
            
            dark = dark
          )
          
        })
      })
    }
  })
  
  output$summary_tbl <- renderTable({
    
    res <- results()
    
    req(res)
    
    res$summary_tbl
    
  }, rownames = FALSE)
}

shinyApp(ui, server)