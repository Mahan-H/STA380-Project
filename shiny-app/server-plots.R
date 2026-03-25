library(ggplot2)
library(tidyr)

model_labels <- c(
  
  noisy = "Noisy QRW",
  
  noiseless = "Noiseless QRW",
  
  classical = "Classical RW"
)

model_cols <- c(
  
  noisy = "#E74C3C",
  
  noiseless = "#2E75B6",
  
  classical = "#27AE60"
)

plot_theme <- function(dark = FALSE) {
  
  if (dark) {
    
    theme_grey(base_size = 13) +
      
      theme(
        
        plot.background  = element_rect(fill = "#212529", colour = NA),
        
        panel.background = element_rect(fill = "#2c3034", colour = NA),
        
        panel.grid.major = element_line(colour = "#495057"),
        
        panel.grid.minor = element_line(colour = "#343a40"),
        
        text             = element_text(colour = "#dee2e6"),
        
        axis.text        = element_text(colour = "#adb5bd"),
        
        axis.ticks       = element_line(colour = "#495057"),
        
        legend.background = element_rect(fill = "#212529", colour = NA),
        
        legend.key       = element_rect(fill = "#2c3034", colour = NA)
      )
    
  } else {
    
    theme_grey(base_size = 13)
  }
}

make_dist_plot_data <- function(dists, selected_models) {
  
  plot_df <- data.frame(
    
    pos = dists$pos,
    
    noisy = dists$noisy,
    
    noiseless = dists$noiseless,
    
    classical = dists$classical
  )
  
  plot_df <- pivot_longer(
    
    plot_df,
    
    cols = c("noisy", "noiseless", "classical"),
    
    names_to = "model",
    
    values_to = "prob"
  )
  
  plot_df <- plot_df[plot_df$model %in% selected_models, , drop = FALSE]
  
  plot_df$model <- factor(
    
    plot_df$model,
    
    levels = selected_models,
    
    labels = model_labels[selected_models]
  )
  
  plot_df
}

make_var_plot_data <- function(var_tbl, selected_models) {
  
  plot_df <- data.frame(
    
    time = var_tbl$time,
    
    noisy = var_tbl$noisy,
    
    noiseless = var_tbl$noiseless,
    
    classical = var_tbl$classical
  )
  
  plot_df <- pivot_longer(
    
    plot_df,
    
    cols = c("noisy", "noiseless", "classical"),
    
    names_to = "model",
    
    values_to = "var"
  )
  
  plot_df <- plot_df[plot_df$model %in% selected_models, , drop = FALSE]
  
  plot_df$model <- factor(
    
    plot_df$model,
    
    levels = selected_models,
    
    labels = model_labels[selected_models]
  )
  
  plot_df
}

plot_dist_overlay <- function(plot_df, selected_models, T, N, channel, p, dark = FALSE) {
  
  sub_text <- paste0(
    
    toupper(substr(channel, 1, 1)),
    
    substr(channel, 2, nchar(channel)),
    
    " channel at p = ", p
  )
  
  p_obj <- plot_ly()
  
  for (m in selected_models) {
    
    lab <- unname(model_labels[m])
    
    col <- unname(model_cols[m])
    
    df_m <- plot_df[plot_df$model == lab, , drop = FALSE]
    
    p_obj <- p_obj |>
      
      add_lines(
        
        data = df_m,
        
        x = ~pos,
        
        y = ~prob,
        
        name = lab,
        
        line = list(color = col, width = 2),
        
        hovertemplate = paste0(
          
          "Model: ", lab,
          
          "<br>Position: %{x}",
          
          "<br>Probability: %{y:.4f}",
          
          "<extra></extra>"
        )
      )
  }
  
  apply_plotly_theme(
    
    p_obj,
    
    title = paste0("Final Position Distributions (T = ", T, ", N = ", N, ")"),
    
    subtitle = sub_text,
    
    x_title = "Position",
    
    y_title = "Probability",
    
    dark = dark,
    
    legend_orientation = "h",
    
    legend_y = 1.12,
    
    show_legend = TRUE
  )
}

plot_var_overlay <- function(plot_df, selected_models, dark = FALSE) {
  
  p_obj <- plot_ly()
  
  for (m in selected_models) {
    
    lab <- unname(model_labels[m])
    
    col <- unname(model_cols[m])
    
    df_m <- plot_df[plot_df$model == lab, , drop = FALSE]
    
    p_obj <- p_obj |>
      
      add_lines(
        
        data = df_m,
        
        x = ~time,
        
        y = ~var,
        
        name = lab,
        
        line = list(color = col, width = 2),
        
        hovertemplate = paste0(
          
          "Model: ", lab,
          
          "<br>Time: %{x}",
          
          "<br>Variance: %{y:.4f}",
          
          "<extra></extra>"
        )
      )
  }
  
  apply_plotly_theme(
    
    p_obj,
    
    title = "Variance vs. Time",
    
    x_title = "Time Step",
    
    y_title = "Var(X_t)",
    
    dark = dark,
    
    legend_orientation = "v",
    
    legend_y = 1,
    
    show_legend = TRUE
  )
}

plot_dist_single <- function(dists, model_name, T, init_coin, dark = FALSE) {
  
  model_titles <- c(
    
    noisy = "Noisy QRW",
    
    noiseless = "Noiseless Hadamard QRW",
    
    classical = "Classical RW"
  )
  
  df_m <- data.frame(
    
    pos = dists$pos,
    
    prob = dists[[model_name]]
  )
  
  subtitle_text <- if (init_coin == "symmetric") {
    
    "Symmetric initial coin state |ψ₀⟩ = (|↑⟩ + i|↓⟩)/√2
    "
  } else {
    
    paste("Initial coin state:", init_coin)
  }
  
  p_obj <- plot_ly(
    
    data = df_m,
    
    x = ~pos,
    
    y = ~prob,
    
    type = "bar",
    
    marker = list(color = unname(model_cols[model_name])),
    
    hovertemplate = paste0(
      
      "Position: %{x}",
      
      "<br>Probability: %{y:.4f}",
      
      "<extra></extra>"
    )
  )
  
  apply_plotly_theme(
    
    p_obj,
    
    title = paste0(model_titles[model_name], " — Final Distribution (T = ", T, ")"),
    
    subtitle = subtitle_text,
    
    x_title = "Position",
    
    y_title = "Probability",
    
    dark = dark,
    
    show_legend = FALSE
  )
}

plot_var_single <- function(var_tbl, model_name, dark = FALSE) {
  
  model_titles <- c(
    
    noisy = "Noisy QRW",
    
    noiseless = "Noiseless Hadamard QRW",
    
    classical = "Classical RW"
  )
  
  df_m <- data.frame(
    
    time = var_tbl$time,
    
    var = var_tbl[[model_name]]
  )
  
  p_obj <- plot_ly(
    
    data = df_m,
    
    x = ~time,
    
    y = ~var,
    
    type = "scatter",
    
    mode = "lines",
    
    line = list(color = unname(model_cols[model_name]), width = 2),
    
    hovertemplate = paste0(
      
      "Time: %{x}",
      
      "<br>Variance: %{y:.4f}",
      
      "<extra></extra>"
    )
  )
  
  apply_plotly_theme(
    
    p_obj,
    
    title = paste0(model_titles[model_name], " — Variance vs. Time"),
    
    x_title = "Time Step",
    
    y_title = "Var(X_t)",
    
    dark = dark,
    
    show_legend = FALSE
  )
}

make_table <- function(results, main_tab, selected_models) {
  
  if (main_tab == "Final Distribution") {
    
    out <- data.frame(pos = results$dists$pos)
    
    if ("noisy" %in% selected_models) out$noisy <- results$dists$noisy
    
    if ("noiseless" %in% selected_models) out$noiseless <- results$dists$noiseless
    
    if ("classical" %in% selected_models) out$classical <- results$dists$classical
    
    return(out)
  }
  
  if (main_tab == "Variance vs Time") {
    
    out <- data.frame(time = results$var_tbl$time)
    
    if ("noisy" %in% selected_models) out$noisy <- results$var_tbl$noisy
    
    if ("noiseless" %in% selected_models) out$noiseless <- results$var_tbl$noiseless
    
    if ("classical" %in% selected_models) out$classical <- results$var_tbl$classical
    
    return(out)
  }
  
  data.frame()
}

