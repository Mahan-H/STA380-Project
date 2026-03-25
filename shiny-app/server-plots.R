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
  
  ggplot(plot_df, aes(x = pos, y = prob, colour = model)) +
    
    geom_line(linewidth = 0.9, alpha = 0.9) +
    
    scale_colour_manual(values = setNames(model_cols[selected_models], model_labels[selected_models])) +
    
    labs(
      
      title = paste0("Final Position Distributions (T = ", T, ", N = ", N, ")"),
      
      subtitle = sub_text,
      
      x = "Position",
      
      y = "Probability",
      
      colour = NULL
    ) +
    
    plot_theme(dark) +
    
    theme(legend.position = "top")
}

plot_var_overlay <- function(plot_df, selected_models, dark = FALSE) {
  
  ggplot(plot_df, aes(x = time, y = var, colour = model)) +
    
    geom_line(linewidth = 0.9) +
    
    scale_colour_manual(values = setNames(model_cols[selected_models], model_labels[selected_models])) +
    
    labs(
      
      title = "Variance vs. Time",
      
      x = "Time Step",
      
      y = "Var(X_t)",
      
      colour = NULL
    ) +
    
    plot_theme(dark) +
    
    theme(legend.position = "right")
}

plot_dist_single <- function(dists, model_name, T, init_coin, dark = FALSE) {
  
  model_titles <- c(
    
    noisy = "Noisy QRW",
    
    noiseless = "Noiseless Hadamard QRW",
    
    classical = "Classical RW"
  )
  
  ggplot(dists, aes(x = pos, y = .data[[model_name]])) +
    
    geom_col(fill = model_cols[model_name], width = 0.8, alpha = 0.85) +
    
    labs(
      
      title = paste0(model_titles[model_name], " — Final Distribution (T = ", T, ")"),
      
      subtitle = if (init_coin == "symmetric") {
        
        "Symmetric initial coin state  |ψ₀⟩ = (|↑⟩ + i|↓⟩)/√2"
        
      } else {
        
        paste("Initial coin state:", init_coin)
      },
      
      x = "Position",
      
      y = "Probability"
    ) +
    
    plot_theme(dark)
}

plot_var_single <- function(var_tbl, model_name, dark = FALSE) {
  
  model_titles <- c(
    
    noisy = "Noisy QRW",
    
    noiseless = "Noiseless Hadamard QRW",
    
    classical = "Classical RW"
  )
  
  ggplot(var_tbl, aes(x = time, y = .data[[model_name]])) +
    
    geom_line(linewidth = 0.9, colour = model_cols[model_name]) +
    
    labs(
      
      title = paste0(model_titles[model_name], " — Variance vs. Time"),
      
      x = "Time Step",
      
      y = "Var(X_t)"
    ) +
    
    plot_theme(dark)
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

