library(testthat)

test_that("pos_grid returns the correct symmetric grid", {
  
  expect_equal(pos_grid(0), 0)
  
  expect_equal(pos_grid(3), -3:3)
  
  expect_length(pos_grid(4), 9)
})

test_that("Hadamard coin has the expected form and is unitary", {
  
  H <- hadamard_coin()
  
  expect_equal(H, (1 / sqrt(2)) * matrix(c(1, 1, 1, -1), nrow = 2, byrow = TRUE))
  
  expect_equal(H %*% Conj(t(H)), diag(2), tolerance = 1e-12)
})

test_that("Pauli matrices square to the identity", {
  
  I2 <- diag(2)
  
  I2_comp <- matrix(c(1+0i, 0+0i, 0+0i, 1+0i), nrow = 2, byrow = TRUE)
  
  expect_equal(pauli_X() %*% pauli_X(), I2, tolerance = 1e-12)
  
  expect_equal(pauli_Y() %*% pauli_Y(), I2_comp, tolerance = 1e-12)
  
  expect_equal(pauli_Z() %*% pauli_Z(), I2, tolerance = 1e-12)
})

test_that("Basis vectors and projectors are correct", {
  
  r <- basis_right()
  
  l <- basis_left()
  
  Pr <- proj_right()
  
  Pl <- proj_left()
  
  expect_equal(r, matrix(c(1, 0), ncol = 1))
  
  expect_equal(l, matrix(c(0, 1), ncol = 1))
  
  expect_equal(Pr, matrix(c(1, 0, 0, 0), nrow = 2, byrow = TRUE))
  
  expect_equal(Pl, matrix(c(0, 0, 0, 1), nrow = 2, byrow = TRUE))
  
  expect_equal(Pr %*% Pr, Pr, tolerance = 1e-12)
  
  expect_equal(Pl %*% Pl, Pl, tolerance = 1e-12)
  
  expect_equal(Pr + Pl, identity_coin(), tolerance = 1e-12)
})

test_that("Initial coin states are constructed correctly", {
  
  expect_equal(init_coin_state("right"), basis_right())
  
  expect_equal(init_coin_state("left"), basis_left())
  
  sym <- init_coin_state("symmetric")
  
  expect_equal(state_norm_sq(sym), 1, tolerance = 1e-12)
  
  expect_equal(sym, (1 / sqrt(2)) * (basis_right() + 1i * basis_left()), tolerance = 1e-12)
})

test_that("initialize_state places the chosen coin state at the origin", {
  
  state <- initialize_state(2, "right")
  
  expect_equal(dim(state), c(2, 5))
  
  expect_equal(state[, 3], c(1, 0))
  
  expect_equal(state_norm_sq(state), 1, tolerance = 1e-12)
  
  state_left <- initialize_state(3, "left")
  
  expect_equal(state_left[, 4], c(0, 1))
})

test_that("state_norm_sq and normalize_state behave correctly", {
  
  state <- matrix(c(2, 0), ncol = 1)
  
  expect_equal(state_norm_sq(state), 4)
  
  state_norm <- normalize_state(state)
  
  expect_equal(state_norm_sq(state_norm), 1, tolerance = 1e-12)
  
  zero_state <- matrix(0, nrow = 2, ncol = 3)
  
  expect_error(normalize_state(zero_state), "nearly 0")
})

test_that("apply_coin_op acts column wise on the coin state", {
  
  state <- initialize_state(1, "right")
  
  out <- apply_coin_op(state, hadamard_coin())
  
  expect_equal(out[, 2], (1 / sqrt(2)) * c(1, 1), tolerance = 1e-12)
})

test_that("apply_shift_op moves right amplitudes right and left amplitudes left", {
  
  state <- matrix(0 + 0i, nrow = 2, ncol = 5)
  
  state[, 3] <- c(1 + 0i, 2 + 0i)
  
  shifted <- apply_shift_op(state)
  
  expect_equal(shifted[1, 4], 1 + 0i)
  
  expect_equal(shifted[2, 2], 2 + 0i)
  
  expect_equal(sum(Mod(shifted)^2), 5, tolerance = 1e-12)
})

test_that("apply_unitary_step produces the expected one-step walk from right state", {
  
  state <- initialize_state(1, "right")
  
  out <- apply_unitary_step(state)
  
  probs <- pos_dist(out)
  
  expect_equal(probs, c(1 / 2, 0, 1 / 2), tolerance = 1e-12)
  
  expect_equal(sum(probs), 1, tolerance = 1e-12)
})

test_that("Kraus operator families satisfy completeness", {
  
  I_comp <- matrix(c(1+0i, 0+0i, 0+0i, 1+0i), nrow = 2, byrow = TRUE)
  
  check_completeness <- function(kraus_ops) {
    
    acc <- matrix(0 + 0i, nrow = 2, ncol = 2)
    
    for (k in seq_along(kraus_ops)) {
      
      acc <- acc + Conj(t(kraus_ops[[k]])) %*% kraus_ops[[k]]
    }
    
    acc
  }
  
  dep <- kraus_dephasing(0.3)
  
  dpol <- kraus_depolarizing(0.6)
  
  expect_equal(check_completeness(dep), I_comp, tolerance = 1e-12)
  
  expect_equal(check_completeness(dpol), I_comp, tolerance = 1e-12)
})

test_that("get_kraus_ops returns the requested family", {
  
  expect_length(get_kraus_ops("dephasing", 0.2), 3)
  
  expect_length(get_kraus_ops("depolarizing", 0.2), 4)
})

test_that("apply_coin_kraus with identity leaves the state unchanged", {
  
  state <- initialize_state(2, "symmetric")
  
  out <- apply_coin_kraus(state, identity_coin())
  
  expect_equal(out, state, tolerance = 1e-12)
})

test_that("kraus_probs are nonnegative and sum to one", {
  
  state <- initialize_state(2, "symmetric")
  
  probs_dep <- kraus_probs(state, kraus_dephasing(0.2))
  
  probs_dpol <- kraus_probs(state, kraus_depolarizing(0.2))
  
  expect_true(all(probs_dep >= 0))
  
  expect_true(all(probs_dpol >= 0))
  
  expect_equal(sum(probs_dep), 1, tolerance = 1e-12)
  
  expect_equal(sum(probs_dpol), 1, tolerance = 1e-12)
})

test_that("sample_kraus_index returns a valid index", {
  
  set.seed(1)
  
  idx <- sample_kraus_index(c(0.1, 0.2, 0.7))
  
  expect_length(idx, 1)
  
  expect_true(idx %in% 1:3)
})

test_that("apply_decoherence_step returns a normalized state of the same dimension", {
  
  set.seed(42)
  
  state <- initialize_state(2, "symmetric")
  
  out <- apply_decoherence_step(state, "dephasing", 0.25)
  
  expect_equal(dim(out), dim(state))
  
  expect_equal(state_norm_sq(out), 1, tolerance = 1e-12)
})

test_that("pos_dist, mean_pos, and var_pos work on simple examples", {
  
  state <- matrix(c(0, 1 / sqrt(2), 0,
                    0, 1i / sqrt(2), 0), nrow = 2, byrow = TRUE)
  
  probs <- pos_dist(state)
  
  expect_equal(probs, c(0, 1, 0), tolerance = 1e-12)
  
  pos <- c(-1, 0, 1)
  
  probs2 <- c(0.25, 0.5, 0.25)
  
  expect_equal(mean_pos(pos, probs2), 0, tolerance = 1e-12)
  
  expect_equal(var_pos(pos, probs2), 0.5, tolerance = 1e-12)
})

test_that("sim_noiseless_qrw returns normalized distributions over time", {
  
  res <- sim_noiseless_qrw(2, "right")
  
  expect_equal(dim(res$dists), c(3, 5))
  
  expect_equal(rowSums(res$dists), c(1, 1, 1), tolerance = 1e-12)
  
  expect_equal(length(res$means), 3)
  
  expect_equal(length(res$vars), 3)
})

test_that("srw_dist_at_t matches small exact classical distributions", {
  
  pos <- pos_grid(2)
  
  expect_equal(srw_dist_at_t(0, pos), c(0, 0, 1, 0, 0), tolerance = 1e-12)
  
  expect_equal(srw_dist_at_t(1, pos), c(0, 0.5, 0, 0.5, 0), tolerance = 1e-12)
  
  expect_equal(srw_dist_at_t(2, pos), c(0.25, 0, 0.5, 0, 0.25), tolerance = 1e-12)
})

test_that("sim_srw returns exact normalized classical distributions", {
  
  res <- sim_srw(3)
  
  expect_equal(dim(res$dists), c(4, 7))
  
  expect_equal(rowSums(res$dists), rep(1, 4), tolerance = 1e-12)
  
  expect_equal(res$vars, 0:3, tolerance = 1e-12)
})

test_that("report table builders return objects of the expected shape", {
  
  noisy <- sim_noisy_qrw(2, 2, "dephasing", "symmetric", 0.1, seed = 1)
  
  noiseless <- sim_noiseless_qrw(2, "symmetric")
  
  classical <- sim_srw(2)
  
  dist_table <- build_dist(noisy, noiseless, classical)
  
  var_table <- build_variance_overlay(noisy, noiseless, classical)
  
  summary_table <- build_summary_table(noisy, noiseless, classical)
  
  expect_equal(names(dist_table), c("pos", "noisy", "noiseless", "classical"))
  
  expect_equal(names(var_table), c("time", "noisy", "noiseless", "classical"))
  
  expect_equal(names(summary_table), c("model", "mean", "variance"))
  
  expect_equal(nrow(summary_table), 3)
})

test_that("sim_noisy_qrw returns averaged normalized distributions and full summaries", {
  
  res <- sim_noisy_qrw(2, 3, "dephasing", "symmetric", 0.1, seed = 1)
  
  expect_equal(rowSums(res$dists), c(1, 1, 1), tolerance = 1e-12)
  
  expect_equal(length(res$means), 3)
  
  expect_equal(length(res$vars), 3)
  
  expect_equal(res$means[1], 0, tolerance = 1e-12)
  
  expect_equal(res$vars[1], 0, tolerance = 1e-12)
})
