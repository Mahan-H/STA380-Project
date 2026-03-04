# STA380 Project 
# Project: Decoherence in a Hadamard Quantum Random Walk on a Line

library(tidyverse)

############################################## CORE QUANTUM OBJECTS / HELPERS

#' Create the position grid for the walk.
#' 
#' @description This returns the integer position grid from \code{-T} to
#' \code{T}, which is used for the quantum and classical walks.
#' 
#' @param T represents the number of walk steps.
#' @return an integer vector of positions from \code{-T} to \code{T}.
#' 
#' @examples
#' pos_grid(2)
#' 
#' @export
pos_grid <- function(T) {
  
  return(seq.int(-T, T))
}

#' Construct the Hadamard coin operator.
#' 
#' @description This returns the \eqn{2 \times 2} Hadamard matrix acting on
#' the coin space. The Hadamard coin is
#' \deqn{
#' H = \frac{1}{\sqrt{2}}
#' \begin{pmatrix}
#' 1 & 1 \\
#' 1 & -1
#' \end{pmatrix}.
#' }
#' 
#' @return a \eqn{2 \times 2} matrix representing the Hadamard coin.
#' 
#' @examples
#' hadamard_coin()
#' 
#' @export
hadamard_coin <- function() {
  
  H <- 1/sqrt(2) * matrix(c(1, 1, 1, -1), nrow = 2, byrow = TRUE)
  
  return(H)
}

#' Construct the Pauli X matrix.
#' 
#' @description This returns the Pauli X operator acting on the
#' two-dimensional coin space.
#' 
#' @return a \eqn{2 \times 2} matrix representing the Pauli X operator.
#' 
#' @examples
#' pauli_X()
#' 
#' @export
pauli_X <- function() {
  
  px <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  
  return(px)
}

#' Construct the Pauli Y matrix.
#' 
#' @description This returns the Pauli Y operator acting on the
#' two-dimensional coin space.
#' 
#' @return a complex \eqn{2 \times 2} matrix representing the Pauli Y operator.
#' 
#' @examples
#' pauli_Y()
#' 
#' @export
pauli_Y <- function() {
  
  py <- matrix(c(0, -1i, 1i, 0), nrow = 2, byrow = TRUE)
  
  return(py)
}

#' Construct the Pauli Z matrix.
#' 
#' @description This returns the Pauli Z operator acting on the
#' two-dimensional coin space.
#' 
#' @return a \eqn{2 \times 2} matrix representing the Pauli Z operator.
#' 
#' @examples
#' pauli_Z()
#' 
#' @export
pauli_Z <- function() {
  
  pz <- matrix(c(1, 0, 0, -1), nrow = 2, byrow = TRUE)
  
  return(pz)
}

#' Construct the right coin basis vector.
#' 
#' @description This returns the right coin basis state as a column vector.
#' 
#' @return a \eqn{2 \times 1} column vector representing the right coin state.
#' 
#' @examples
#' basis_right()
#' 
#' @export
basis_right <- function() {
  
  matrix(c(1, 0), ncol = 1)
}

#' Construct the left coin basis vector.
#' 
#' @description This returns the left coin basis state as a column vector.
#' 
#' @return a \eqn{2 \times 1} column vector representing the left coin state.
#' 
#' @examples
#' basis_left()
#' 
#' @export
basis_left <- function() {
  
  matrix(c(0, 1), ncol = 1)
}

#' Construct the projector onto the right coin state.
#' 
#' @description This returns the rank one projector
#' \eqn{|\mathrm{right}\rangle\langle\mathrm{right}|}.
#' 
#' @return a \eqn{2 \times 2} projector matrix.
#' 
#' @examples
#' proj_right()
#' 
#' @export
proj_right <- function() {
  
  # |right>
  r <- basis_right()
  
  # |right><right|
  proj <- r %*% Conj(t(r))
  
  return(proj)
}

#' Construct the projector onto the left coin state.
#' 
#' @description This returns the rank one projector
#' \eqn{|\mathrm{left}\rangle\langle\mathrm{left}|}.
#' 
#' @return a \eqn{2 \times 2} projector matrix.
#' 
#' @examples
#' proj_left()
#' 
#' @export
proj_left <- function() {
  
  # |left>
  l <- basis_left()
  
  # |left><left|
  proj <- l %*% Conj(t(l))
  
  return(proj)
}

#' Construct the identity operator on the coin space.
#' 
#' @description This returns the identity matrix acting on the
#' two-dimensional coin space.
#' 
#' @return a \eqn{2 \times 2} identity matrix.
#' 
#' @examples
#' identity_coin()
#' 
#' @export
identity_coin <- function() {
  
  I <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  
  return(I)
}

#################################################### INITIAL STATE CONSTRUCTION

#' Construct the initial coin state.
#' 
#' @description This returns the chosen initial coin state for the walk.
#' The supported choices are \code{"right"}, \code{"left"}, and
#' \code{"symmetric"}.
#' 
#' @param init_coin represents the initial coin state used to start the walk.
#' @return a \eqn{2 \times 1} column vector representing the initial coin state.
#' 
#' @examples
#' init_coin_state("right")
#' init_coin_state("symmetric")
#' 
#' @export
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

#' Initialize the full quantum walk state.
#' 
#' @description This creates the state matrix for the walk, placing the walker
#' at the origin with the chosen initial coin state.
#' 
#' @param T represents the number of walk steps.
#' @param init_coin represents the initial coin state used to start the walk.
#' @return a \eqn{2 \times (2T+1)} matrix representing the initial walk state.
#' 
#' @examples
#' initialize_state(2, "right")
#' 
#' @export
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

#' Compute the squared norm of a walk state.
#' 
#' @description This computes the total probability mass of the walk state by
#' summing the squared moduli of all amplitudes.
#' 
#' @param state represents the state matrix of the walk.
#' @return a numeric value equal to the squared norm of the state.
#' 
#' @examples
#' state <- initialize_state(1, "right")
#' state_norm_sq(state)
#' 
#' @export
state_norm_sq <- function(state) {
  
  # |amplitude|^2
  norm <- sum(Mod(state)^2)
  
  return(norm)
}

#' Normalize a walk state.
#' 
#' @description This rescales the walk state so that its squared norm is 1.
#' 
#' @param state represents the state matrix of the walk.
#' @param tol represents the tolerance used to detect a near-zero norm.
#' @return a normalized state matrix.
#' 
#' @examples
#' state <- initialize_state(1, "right")
#' normalize_state(state)
#' 
#' @export
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

#' Apply a coin operator to a walk state.
#' 
#' @description This applies a \eqn{2 \times 2} coin operator independently
#' at each position of the walk.
#' 
#' @param state represents the state matrix of the walk.
#' @param coin_matrix a \eqn{2 \times 2} matrix acting on the coin space.
#' @return the updated state matrix after the coin operator is applied.
#' 
#' @examples
#' state <- initialize_state(1, "right")
#' apply_coin_op(state, hadamard_coin())
#' 
#' @export
apply_coin_op <- function(state, coin_matrix) {
  
  for(i in 1:ncol(state)){
    
    state[, i] <- coin_matrix %*% state[, i]
  }
  
  return(state)
}

#' Apply the conditional shift operator.
#' 
#' @description This applies the conditional shift operator for the walk. 
#' The shift operator sends the right coin state one step to the
#' right and the left coin state one step to the left:
#' \deqn{
#' S|x,\mathrm{right}\rangle = |x+1,\mathrm{right}\rangle,
#' \qquad
#' S|x,\mathrm{left}\rangle = |x-1,\mathrm{left}\rangle.
#' }
#' 
#' @param state represents the state matrix of the walk.
#' @return the updated state matrix after the conditional shift.
#' 
#' @examples
#' state <- initialize_state(1, "right")
#' apply_shift_op(state)
#' 
#' @export
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

#' Apply one noiseless Hadamard quantum walk step.
#' 
#' @description This applies one noiseless step of the Hadamard quantum random
#' walk. The unitary step is
#' \deqn{
#' U = S(H \otimes I_p),
#' }
#' where \eqn{H} is the Hadamard coin, \eqn{S} is the conditional shift, and
#' \eqn{I_p} is the identity on the position space.
#' 
#' @param state represents the state matrix of the walk.
#' @return the updated state matrix after one noiseless walk step.
#' 
#' @examples
#' state <- initialize_state(2, "symmetric")
#' apply_unitary_step(state)
#' 
#' @export
apply_unitary_step <- function(state) {
  
  H <- hadamard_coin()
  
  state <- apply_coin_op(state, H)
  
  new_state <- apply_shift_op(state)
  
  return(new_state)
}

######################################### KRAUS OPERATORS / DECOHERENCE CHANNELS

#' Construct the Kraus operators for the dephasing channel.
#' 
#' @description This returns the Kraus operator family for the dephasing
#' channel acting on the coin space. The operators are
#' \deqn{
#' K_0 = \sqrt{1-p}\,I_c,\qquad
#' K_1 = \sqrt{p}\,|\mathrm{right}\rangle\langle\mathrm{right}|,\qquad
#' K_2 = \sqrt{p}\,|\mathrm{left}\rangle\langle\mathrm{left}|.
#' }
#' 
#' @param p represents the decoherence strength.
#' @return a list of Kraus operators for the dephasing channel.
#' 
#' @examples
#' kraus_dephasing(0.1)
#' 
#' @export
kraus_dephasing <- function(p) {
  
  k0 <- sqrt(1 - p) * identity_coin()
  
  k1 <- sqrt(p) * proj_right()
  
  k2 <- sqrt(p) * proj_left()
  
  kraus_ops <- list(k0, k1, k2)
  
  return(kraus_ops)
}

#' Construct the Kraus operators for the depolarizing channel.
#' 
#' @description This returns the Kraus operator family for the depolarizing
#' channel acting on the coin space. The operators are
#' \deqn{
#' K_0 = \sqrt{1-p}\,I_c,\qquad
#' K_1 = \sqrt{p/3}\,X,\qquad
#' K_2 = \sqrt{p/3}\,Y,\qquad
#' K_3 = \sqrt{p/3}\,Z.
#' }
#' 
#' @param p represents the decoherence strength.
#' @return a list of Kraus operators for the depolarizing channel.
#' 
#' @examples
#' kraus_depolarizing(0.1)
#' 
#' @export
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

#' Select the Kraus operator family for a channel.
#' 
#' @description This returns the appropriate Kraus operator family for either
#' the dephasing or depolarizing channel.
#' 
#' @param channel represents the selected decoherence channel.
#' @param p represents the decoherence strength.
#' @return a list of Kraus operators for the selected channel.
#' 
#' @examples
#' get_kraus_ops("dephasing", 0.1)
#' 
#' @export
get_kraus_ops <- function(channel, p) {
  
  if(channel == "dephasing"){
    
    return(kraus_dephasing(p))
    
  } else if(channel == "depolarizing") {
    
    return(kraus_depolarizing(p))
    
  } else {
    
    stop("Invalid channel.")
  }
}

#' Apply a Kraus operator to the coin state.
#' 
#' @description This applies a single Kraus operator to the coin component
#' at every position of the walk. The full operator on the joint space is
#' \deqn{
#' \widetilde{K} = K \otimes I_p,
#' }
#' where \eqn{K} acts on the coin space and \eqn{I_p} acts on the position
#' space.
#' 
#' @param state represents the state matrix of the walk.
#' @param K a \eqn{2 \times 2} Kraus operator acting on the coin space.
#' @return the updated state matrix after the Kraus operator is applied.
#' 
#' @examples
#' state <- initialize_state(1, "right")
#' ks <- kraus_dephasing(0.1)
#' apply_coin_kraus(state, ks[[1]])
#' 
#' @export
apply_coin_kraus <- function(state, K) {
  
  for(i in 1:ncol(state)){
    
    state[, i] <- K %*% state[, i]
  }
  
  return(state)
}

#' Compute Kraus outcome probabilities.
#' 
#' @description This computes the probabilities associated with each Kraus
#' operator for the current walk state. The probability of outcome \eqn{k} is
#' \deqn{
#' p_k = \|\widetilde{K}_k |\phi\rangle\|^2,
#' }
#' where \eqn{|\phi\rangle} is the state after the unitary step.
#' 
#' @param state represents the state matrix of the walk.
#' @param kraus_ops a list of Kraus operators.
#' @param tol represents the tolerance used for small numerical values.
#' @return a numeric vector of Kraus outcome probabilities.
#' 
#' @examples
#' state <- initialize_state(1, "symmetric")
#' ks <- kraus_dephasing(0.1)
#' kraus_probs(state, ks)
#' 
#' @export
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

#' Sample a Kraus index.
#' 
#' @description This samples one Kraus operator index according to the
#' supplied probability vector.
#' 
#' @param probs a numeric probability vector.
#' @return an integer giving the sampled Kraus index.
#' 
#' @examples
#' set.seed(1)
#' sample_kraus_index(c(0.2, 0.8))
#' 
#' @export
sample_kraus_index <- function(probs) {
  
  samp_index <- sample(seq(1:length(probs)), size = 1, prob = probs)
  
  return(samp_index)
}

#' Apply one decoherence step.
#' 
#' @description This applies one sampled decoherence step to the current walk
#' state and then normalizes the result. After sampling outcome \eqn{k}, the updated state is
#' \deqn{
#' |\psi\rangle \leftarrow
#' \frac{\widetilde{K}_k |\phi\rangle}{\sqrt{p_k}},
#' }
#' where \eqn{p_k = \|\widetilde{K}_k |\phi\rangle\|^2}.
#' 
#' @param state represents the state matrix of the walk.
#' @param channel represents the selected decoherence channel.
#' @param p represents the decoherence strength.
#' @return the updated normalized state matrix after one decoherence step.
#' 
#' @examples
#' set.seed(1)
#' state <- initialize_state(1, "symmetric")
#' apply_decoherence_step(state, "dephasing", 0.1)
#' 
#' @export
apply_decoherence_step <- function(state, channel, p) {
  
  kraus_ops <- get_kraus_ops(channel, p)
  
  probs <- kraus_probs(state, kraus_ops)
  
  samp_index <- sample_kraus_index(probs)
  
  new_state <- apply_coin_kraus(state, kraus_ops[[samp_index]])
  
  normalized_state <- normalize_state(new_state)
  
  return(normalized_state)
}

########################################## MONTE CARLO QUANTUM TRAJECTORIES


#' Simulate a noisy Hadamard quantum random walk.
#' 
#' @description This simulates a decohering Hadamard quantum random walk on a
#' line using Monte Carlo trajectories. The Monte Carlo estimate of the position distribution is
#' \deqn{
#' \widehat{P}_t(x) = \frac{1}{N}\sum_{j=1}^N P_t^{(j)}(x),
#' }
#' where \eqn{N} is the number of trajectories.
#' 
#' @param T represents the number of walk steps.
#' @param N represents the number of Monte Carlo trajectories.
#' @param channel represents the selected decoherence channel. The allowed
#' choices are \code{"dephasing"} and \code{"depolarizing"}.
#' @param init_coin represents the initial coin state. The allowed choices are
#' \code{"right"}, \code{"left"}, and \code{"symmetric"}.
#' @param p represents the decoherence strength.
#' @param seed an optional random seed used for reproducibility.
#' @return a list containing the number of steps, the time indexed position
#' distributions, the mean positions, and the variances.
#' 
#' @examples
#' set.seed(1)
#' sim_noisy_qrw(T = 3, N = 5, channel = "dephasing",
#'               init_coin = "symmetric", p = 0.1, seed = 1)
#'               
#' @export
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
  
  normalize_dists <- dists/N
  
  means <- numeric(length(times))
  
  vars <- numeric(length(times))
  
  for(t in 0:T){
    
    probs <- normalize_dists[t + 1, ]
    
    means[t + 1] <- mean_pos(pos, probs)
    
    vars[t + 1] <- var_pos(pos, probs)
  }
  
  data <- list(T = T, dists = normalize_dists, means = means, vars = vars)
  
  return(data)
}

################################ POSITION DISTRIBUTION

#' Compute the position distribution of a walk state.
#' 
#' @description This computes the position probabilities by summing the right
#' and left coin probabilities at each position. For one trajectory the position distribution is
#' \deqn{
#' P_t^{(j)}(x) = |a_t^{(j)}(x)|^2 + |b_t^{(j)}(x)|^2.
#' }
#' 
#' @param state represents the state matrix of the walk.
#' @return a numeric vector of position probabilities.
#' 
#' @examples
#' state <- initialize_state(2, "right")
#' pos_dist(state)
#' 
#' @export
pos_dist <- function(state) {
  
  n <- ncol(state)
  
  probs <- numeric(n)
  
  for(i in 1:n){
    
    prob <- abs(state[1, i])^2 + abs(state[2, i])^2
    
    probs[i] <- prob
  }
  
  return(probs)
}

#' Compute the mean position.
#' 
#' @description This computes the mean of a probability distribution on
#' the position grid. The mean is
#' \deqn{
#' \mu_t = \sum_x x\,P_t(x).
#' }
#' 
#' @param pos represents the position grid.
#' @param probs represents the probability distribution on the position grid.
#' @return a numeric value representing the mean position.
#' 
#' @examples
#' mean_pos(c(-1, 0, 1), c(0.25, 0.5, 0.25))
#' 
#' @export
mean_pos <- function(pos, probs) {
  
  mean <- sum(pos * probs)
  
  return(mean)
}

#' Compute the position variance.
#' 
#' @description This computes the variance of a probability distribution on
#' the position grid. The variance is
#' \deqn{
#' \mathrm{Var}(X_t) = \sum_x (x - \mu_t)^2 P_t(x).
#' }
#' 
#' @param pos represents the position grid.
#' @param probs represents the probability distribution on the position grid.
#' @return a numeric value representing the position variance.
#' 
#' @examples
#' var_pos(c(-1, 0, 1), c(0.25, 0.5, 0.25))
#' 
#' @export
var_pos <- function(pos, probs) {
  
  mean <- mean_pos(pos, probs)
  
  var <- sum((pos - mean)^2 * probs)
  
  return(var)
}


######################################################## BASELINE

#' Simulate a noiseless Hadamard quantum random walk.
#' 
#' @description This simulates the deterministic Hadamard quantum random walk
#' on a line. The evolution is generated by repeated application of
#' \deqn{
#' U = S(H \otimes I_p).
#' }
#' 
#' @param T represents the number of walk steps.
#' @param init_coin represents the initial coin state. The allowed choices are
#' \code{"right"}, \code{"left"}, and \code{"symmetric"}.
#' @return a list containing the number of steps, the time indexed position
#' distributions, the mean positions, and the variances.
#' 
#' @examples
#' sim_noiseless_qrw(T = 3, init_coin = "symmetric")
#' 
#' @export
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

#' Compute the symmetric random walk distribution at time t.
#' 
#' @description This computes the exact distribution of a symmetric classical
#' random walk after \code{t} steps. If the walker is at position \eqn{x} after
#' \eqn{t} steps, then the number of right steps is
#' \deqn{
#' k = \frac{t + x}{2},
#' }
#' and the probability is
#' \deqn{
#' P(X_t = x) = \binom{t}{k}\left(\frac{1}{2}\right)^t,
#' }
#' whenever \eqn{|x| \le t} and \eqn{t + x} is even.
#' 
#' @param t represents the current time step.
#' @param pos represents the position grid.
#' @return a numeric vector of probabilities on the position grid.
#' 
#' @examples
#' srw_dist_at_t(2, pos_grid(2))
#' 
#' @export
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

#' Simulate a symmetric classical random walk.
#' 
#' @description This computes the exact position distributions, mean positions,
#' and variances for a symmetric classical random walk on a line from time
#' \code{0} to time \code{T}. For the symmetric random walk,
#' \eqn{\mathbb{E}[X_t] = 0} and \eqn{\mathrm{Var}(X_t) = t}.
#' 
#' @param T represents the number of walk steps.
#' @return a list containing the number of steps, the time indexed position
#' distributions, the mean positions, and the variances.
#' 
#' @examples
#' sim_srw(3)
#' 
#' @export
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

#' Build the final distribution comparison table.
#' 
#' @description This combines the final-time distributions of the noisy quantum
#' walk, the noiseless quantum walk, and the classical random walk into one
#' data frame.
#' 
#' @param noisy_w represents the output of \code{sim_noisy_qrw()}.
#' @param noiseless_w represents the output of \code{sim_noiseless_qrw()}.
#' @param srw represents the output of \code{sim_srw()}.
#' @return a data frame containing the position grid and the final time
#' distributions for the three models.
#' 
#' @examples
#' noisy <- sim_noisy_qrw(2, 5, "dephasing", "symmetric", 0.1, seed = 1)
#' noiseless <- sim_noiseless_qrw(2, "symmetric")
#' classical <- sim_srw(2)
#' build_dist(noisy, noiseless, classical)
#' 
#' @export
build_dist <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_w$T
  
  pos <- pos_grid(T)
  
  noisy_dist <- noisy_w$dists[T + 1,]
  
  noiseless_dist <- noiseless_w$dists[T + 1,]
  
  srw_dist <- srw$dists[T + 1,]
  
  model <- data.frame(pos = pos, noisy_dist, noiseless_dist, srw_dist)
  
  names(model)[2] <- "noisy"
  
  names(model)[3] <- "noiseless"
  
  names(model)[4] <- "classical"
  
  return(model)
}

#' Build the variance comparison table.
#' 
#' @description This combines the variance over time outputs of the noisy
#' quantum walk, the noiseless quantum walk, and the classical random walk
#' into one data frame.
#' 
#' @param noisy_w represents the output of \code{sim_noisy_qrw()}.
#' @param noiseless_w represents the output of \code{sim_noiseless_qrw()}.
#' @param srw represents the output of \code{sim_srw()}.
#' @return a data frame containing the time grid and the variance values for
#' the three models.
#' 
#' @examples
#' noisy <- sim_noisy_qrw(2, 5, "dephasing", "symmetric", 0.1, seed = 1)
#' noiseless <- sim_noiseless_qrw(2, "symmetric")
#' classical <- sim_srw(2)
#' build_variance_overlay(noisy, noiseless, classical)
#' 
#' @export
build_variance_overlay <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_w$T
  
  time <- 0:T
  
  noisy_vars <- noisy_w$vars
  
  noiseless_vars <- noiseless_w$vars
  
  srw_vars <- srw$vars
  
  model <- data.frame(time = time, noisy_vars, noiseless_vars, srw_vars)
  
  names(model)[2] <- "noisy"
  
  names(model)[3] <- "noiseless"
  
  names(model)[4] <- "classical"
  
  return(model)
}

#' Build a summary table comparing walk models
#'
#' Creates a final time summary table comparing the noisy quantum random
#' walk, noiseless quantum random walk, and symmetric classical random walk.
#'
#' @param noisy_w Result object returned by \code{sim_noisy_qrw()}.
#' @param noiseless_w Result object returned by \code{sim_noiseless_qrw()}.
#' @param srw Result object returned by \code{sim_srw()}.
#'
#' @return A data frame with one row per model and columns for final-time
#' mean and variance.
#'
#' @examples
#' noisy <- sim_noisy_qrw(4, 20, "dephasing", "symmetric", 0.1, seed = 1)
#' noiseless <- sim_noiseless_qrw(4, "symmetric")
#' classical <- sim_srw(4)
#' build_summary_table(noisy, noiseless, classical)
#'
#' @export
build_summary_table <- function(noisy_w, noiseless_w, srw) {
  
  T <- noisy_w$T
  
  noisy_mean <- noisy_w$means[T + 1]
  
  noisy_var <- noisy_w$vars[T + 1]
  
  noiseless_mean <- noiseless_w$means[T + 1]
  
  noiseless_var <- noiseless_w$vars[T + 1]
  
  srw_mean <- srw$means[T + 1]
  
  srw_var <- srw$vars[T + 1]
  
  summary <- data.frame(c("noisy", "noiseless", "classical"),
                        c(noisy_mean, noiseless_mean, srw_mean),
                        c(noisy_var, noiseless_var, srw_var))
  
  names(summary)[1] <- "model"
  
  names(summary)[2] <- "mean"
  
  names(summary)[3] <- "variance"
  
  return(summary)
}
