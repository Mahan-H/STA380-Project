# STA380 Project 
# Project: Decoherence in a Hadamard Quantum Random Walk on a Line

############################################## CORE QUANTUM OBJECTS / HELPERS

pos_grid <- function(T) {
  
  return(seq.int(-T, T))
}

hadamard_coin <- function() {
  
  H <- 1/sqrt(2) * matrix(c(1, 1, 1, -1), nrow = 2, byrow = TRUE)
  
  return(H)
}

pauli_X <- function() {
  
  px <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  
  return(px)
}

pauli_Y <- function() {
  
  py <- matrix(c(0, -1i, 1i, 0), nrow = 2, byrow = TRUE)
  
  return(py)
}

pauli_Z <- function() {
  
  pz <- matrix(c(1, 0, 0, -1), nrow = 2, byrow = TRUE)
  
  return(pz)
}

basis_right <- function() {
  matrix(c(1, 0), ncol = 1)
}

basis_left <- function() {
  matrix(c(0, 1), ncol = 1)
}

proj_right <- function() {
  
  # |right>
  r <- basis_right()
  
  # |right><right|
  proj <- r %*% Conj(t(r))
  
  return(proj)
}

proj_left <- function() {
  
  # |left>
  l <- basis_left()
  
  # |left><left|
  proj <- l %*% Conj(t(l))
  
  return(proj)
}

identity_coin <- function() {
  
  I <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  
  return(I)
}

#################################################### INITIAL STATE CONSTRUCTION

init_coin_state <- function(init_coin) {
  
  if(init_coin == "right"){
    
    # |right> = (1, 0)^T
    return(basis_right())
    
  } else if (init_coin == "left"){
    
    # |left> = (0, 1)^T
    return(basis_left())
    
  } else if(init_coin == "symmetric"){
    
    # symmetric = 1/sqrt(2) * (|right> + i * |left>)
    sym <- 1/sqrt(2) * (basis_right() + 1i * basis_left())
    
    return(sym)
  }
}

initialize_state <- function(T, init_coin) {
  
  # -T,...,0,...,T
  pos <- pos_grid(T)
  
  # rows are the coin states, columns are the positions
  state <- matrix(0, nrow = 2, ncol = length(pos))
  
  # initial coin state
  c0 <- init_coin_state(init_coin)
  
  # starting position on the grid
  mid <- which(pos == 0)
  
  state[, mid] <- c0
  
  return(state)
}


state_norm_sq <- function(state) {
  
  # |amplitude|^2
  norm <- sum(Mod(state)^2)
  
  return(norm)
}

normalize_state <- function(state, tol = 1e-12) {
  
  norm <- state_norm_sq(state)
  
  # if the norm is nearly 0, then it will cause a division by 0 error during the normalization
  if (norm < tol) {
    
    stop("Cannot normalize since the norm is nearly 0.")
    
  } else {
    
    # normalize the state
    state <- state / sqrt(norm)
  }
  
  return(state)
}

############################################# NOISELESS HADAMARD WALK OPERATORS

apply_coin_op <- function(state, coin_matrix) {
  
  for(i in 1:ncol(state)){
    
    state[, i] <- coin_matrix %*% state[, i]
  }
  
  return(state)
}


apply_shift_op <- function(state) {
  
  n <- ncol(state)
  
  new_state <- matrix(0, nrow = 2, ncol = n)
  
  for(i in 1:n){
    
    if(i < n){
      
      new_state[1, i + 1] <- state[1, i]
    }
    
    if(i - 1 > 0){
      
      new_state[2, i - 1] <- state[2, i]
    }
  }
  
  return(new_state)
}

# U=S(H ⊗ I_p)
apply_unitary_step <- function(state) {
  
  H <- hadamard_coin()
  
  state <- apply_coin_op(state, H)
  
  new_state <- apply_shift_op(state)
  
  return(new_state)
}

######################################### KRAUS OPERATORS / DECOHERENCE CHANNELS

kraus_dephasing <- function(p) {
  
  k0 <- sqrt(1 - p) * identity_coin()
  
  k1 <- sqrt(p) * proj_right()
  
  k2 <- sqrt(p) * proj_left()
  
  kraus_ops <- list(k0, k1, k2)
  
  return(kraus_ops)
}

kraus_depolarizing <- function(p) {
  
  x <- pauli_X()
  
  y <- pauli_Y()
  
  z <- pauli_Z()
  
  k0 <- sqrt(1 - p) * identity_coin()
  
  k1 <- sqrt(p/3) * x
  
  k2 <- sqrt(p/3) * y
  
  k3 <- sqrt(p/3) * z
  
  kraus_ops <- list(k0, k1, k2, k3)
  
  return(kraus_ops)
}

get_kraus_ops <- function(channel, p) {
  
  if(channel == "dephasing"){
    
    return(kraus_dephasing(p))
    
  } else if(channel == "depolarizing") {
    
    return(kraus_depolarizing(p))
    
  } else {
    
    stop("Invalid channel.")
  }
}


# K˜ = K ⊗ I_p
apply_coin_kraus <- function(state, K) {
  
  for(i in 1:ncol(state)){
    
    state[, i] <- K %*% state[, i]
  }
  
  return(state)
}

# ∥K˜_k|ϕ⟩∥^2
kraus_probs <- function(state, kraus_ops, tol = 1e-12) {
  
  n <- length(kraus_ops)
  
  probs <- numeric(n)
  
  for(k in 1:n){
    
    new_state <- apply_coin_kraus(state, kraus_ops[[k]])
    
    pk <- state_norm_sq(new_state)
    
    probs[k] <- pk
  }
  
  # remove tiny negative values caused by numerical roundoff, could cause errors if left there since the probabilities must be non negative
  probs[abs(probs) < tol] <- 0
  
  total_prob <- sum(probs)
  
  # if the total probability is nearly 0, then it will cause a division by 0 error during the normalization
  if (total_prob < tol) {
    stop("Total Kraus probability is neary 0.")
  }
  
  # renormalize in case of small numerical error
  probs <- probs / total_prob
  
  return(probs)
}

sample_kraus_index <- function(probs) {
  
  samp_index <- sample(seq(1:length(probs)), size = 1, prob = probs)
  
  return(samp_index)
}

apply_decoherence_step <- function(state, channel, p) {
  
  kraus_ops <- get_kraus_ops(channel, p)
  
  probs <- kraus_probs(state, kraus_ops)
  
  samp_index <- sample_kraus_index(probs)
  
  new_state <- apply_coin_kraus(state, kraus_ops[[samp_index]])
  
  normalized_state <- normalize_state(new_state)
  
  return(normalized_state)
}

########################################## MONTE CARLO QUANTUM TRAJECTORIES

sim_noisy_qrw <- function(T, N, channel, init_coin, p, seed = NULL) {
  
  if(!is.null(seed)){
    
    set.seed(seed)
  }
  
  pos <- pos_grid(T)
  
  times <- 0:T
  
  dists <- matrix(0, nrow = length(times), ncol = length(pos))
  
  for(i in 1:N){
    
    state <- initialize_state(T, init_coin)
    
    dists[1, ] <- dists[1, ] + pos_dist(state)
    
    for(t in 1:T){
      
      state <- apply_unitary_step(state)
      
      state <- apply_decoherence_step(state, channel, p)
      
      dists[t + 1, ] <- dists[t + 1, ] + pos_dist(state)
    }
  }
  
  avg_dists <- dists/N
  
  means <- numeric(length(times))
  
  vars <- numeric(length(times))
  
  for(t in 1:T){
    
    probs <- avg_dists[t, ]
    
    means[t] <- mean_pos(pos, probs)
    
    vars[t] <- var_pos(pos, probs)
  }
  
  data <- list(T = T, dists = dists, means = means, vars = vars)
  
  return(data)
}

################################ POSITION DISTRIBUTION

# P_t^(j)(x) = |a_t^(j)(x)|^2 + |b_t^(j)(x)|^2
pos_dist <- function(state) {
  
  n <- ncol(state)
  
  probs <- numeric(n)
  
  for(i in 1:n){
    
    prob <- abs(state[1, i])^2 + abs(state[2, i])^2
    
    probs[i] <- prob
  }
  
  return(probs)
}

mean_pos <- function(pos, probs) {
  
  mean <- sum(pos * probs)
  
  return(mean)
}

var_pos <- function(pos, probs) {
  
  mean <- mean_pos(pos, probs)
  
  var <- sum((pos - mean)^2 * probs)
  
  return(var)
}


######################################################## BASELINE


sim_noiseless_qrw <- function(T, init_coin) {
  
  pos <- pos_grid(T)
  
  times <- 0:T
  
  dists <- matrix(0, nrow = length(times), ncol = length(pos))
  
  means <- numeric(length(times))
  
  vars <- numeric(length(times))
  
  state <- initialize_state(T, init_coin)
  
  p0 <- pos_dist(state)
  
  dists[1, ] <- p0
  
  means[1] <- mean_pos(pos, p0)
  
  vars[1] <- var_pos(pos, p0)
  
  for(t in 1:T){
    
    state <- apply_unitary_step(state)
    
    prob <- pos_dist(state)
    
    dists[t + 1, ] <- prob
    
    means[t + 1] <- mean_pos(pos, prob)
    
    vars[t + 1] <- var_pos(pos, prob)
  }
  
  data <- list(T = T, dists = dists, means = means, vars = vars)
  
  return(data)
}

srw_dist_at_t <- function(t, pos) {
  
  n <- length(pos)
  
  probs <- numeric(n)
  
  for(i in 1:n){
    
    if(abs(pos[i]) <= t & (t + pos[i]) %% 2 == 0){
      
      k = (t + pos[i])/2
      
      probs[i] <-  choose(t, k) * (1/2)^t
    }
  }
  
  return(probs)
}

sim_srw <- function(T) {
  
  pos <- pos_grid(T)
  
  times <- 0:T
  
  dists <- matrix(0, nrow = length(times), ncol = length(pos))
  
  means <- numeric(length(times))
  
  vars <- numeric(length(times))
  
  for(t in 0:T){
    
    prob <- srw_dist_at_t(t, pos)
    
    dists[t + 1, ] <- prob
    
    means[t + 1] <- mean_pos(pos, prob)
    
    vars[t + 1] <- var_pos(pos, prob)
  }
  
  data <- list(T = T, dists = dists, means = means, vars = vars)
  
  return(data)
}

################################################### REPORT TABLES

build_dist <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_dist$T
  
  pos <- pos_grid(T)
  
  noisy <- data.frame(pos = pos, noisy_w$dists[t + 1,], model = "noisy")
  
  noiseless <- data.frame(pos = pos, noiseless_w$dists[t + 1,], model = "noiseless")
  
  srw <- data.frame(pos = pos, srw$dists[t + 1,], model = "classical")
  
  models <- rbind(noisy, noiseless, srw)
  
  return(models)
}

build_variance_overlay <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_dist$T
  
  time <- 0:T
  
  noisy <- data.frame(time = time, noisy_w$vars, model = "noisy")
  
  noiseless <- data.frame(time = time, noiseless_w$vars, model = "noiseless")
  
  srw <- data.frame(time = time, srw$vars, model = "classical")
  
  models <- rbind(noisy, noiseless, srw)
  
  return(models)
}

build_summary_table <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_w$T
  
  noisy_mean <- noisy_w$means[T + 1]
  
  noisy_var <- noisy_w$vars[T + 1]
  
  noiseless_mean <- noiseless_w$means[T + 1]
  
  noiseless_var <- noiseless_w$vars[T + 1]
  
  srw_mean <- srw$means[T + 1]
  
  srw_var <- srw$vars[T + 1]
  
  summary <- data.frame(c(noisy_mean, noiseless_mean, srw_mean), c(noisy_var, noiseless_var, srw_var))
}



############################ PLOTTING OUTPUT


plot_dist_overlay <- function(dist_overlay_df) {
  # TODO:
}

plot_var_overlay <- function(var_overlay_df) {
  # TODO:
}