# Decoherence in a Hadamard Quantum Random Walk on a Line

An R-based simulation project for studying how **decoherence changes the behaviour of a one-dimensional Hadamard quantum random walk**.

The project compares three models:

- **Noiseless Hadamard quantum random walk**
- **Noisy Hadamard quantum random walk** with dephasing or depolarizing noise
- **Classical symmetric random walk**

The noisy quantum walk is simulated using **Monte Carlo quantum trajectories** and Kraus operators. An interactive **Shiny application** is included for exploring how the walk changes as the decoherence strength, noise channel, number of steps, and initial coin state are varied.

## Live Shiny App

Explore the simulation interactively:

**https://019cf4b9-90ba-d132-60fa-085607aa91d5.share.connect.posit.cloud**

---

## Project Overview

A classical random walk evolves through probabilistic left/right moves. A quantum random walk instead evolves through amplitudes, superposition, and interference.

For a coherent Hadamard quantum walk, these interference effects generally produce a probability distribution that spreads much faster than a classical random walk. In particular, the variance of an ideal quantum walk typically grows ballistically,

$$
\operatorname{Var}(X_t) \propto t^2,
$$

while a classical symmetric random walk has diffusive spreading,

$$
\operatorname{Var}(X_t) = t.
$$

Real quantum systems are not perfectly isolated. Interaction with an environment destroys phase information and produces **decoherence**. This project investigates that transition by introducing noise into the coin state and comparing the resulting walk with both the noiseless quantum model and the classical benchmark.

---

## Mathematical Model

### 1. Hilbert space

The quantum walker has two components:

1. a **coin state**
2. a **position state**

The coin basis is

$$
|\mathrm{right}\rangle =
\begin{pmatrix}
1\\
0
\end{pmatrix},
\qquad
|\mathrm{left}\rangle =
\begin{pmatrix}
0\\
1
\end{pmatrix}.
$$

The total state belongs to the tensor-product space

$$
\mathcal{H}
=
\mathcal{H}_{\text{coin}}
\otimes
\mathcal{H}_{\text{position}}.
$$

---

### 2. Hadamard coin

At each noiseless quantum-walk step, the Hadamard operator acts on the coin:

$$
H =
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1 & 1\\
1 & -1
\end{pmatrix}.
$$

The Hadamard coin creates a superposition of left- and right-moving amplitudes.

---

### 3. Conditional shift

After the coin operation, the shift operator moves the walker according to the coin state:

$$
S|x,\mathrm{right}\rangle
=
|x+1,\mathrm{right}\rangle,
$$

$$
S|x,\mathrm{left}\rangle
=
|x-1,\mathrm{left}\rangle.
$$

One noiseless quantum-walk step is therefore

$$
U = S(H \otimes I_p),
$$

where \(I_p\) is the identity operator on position space.

---

## Initial Coin States

The implementation supports three initial coin states:

### Right

$$
|\psi_c\rangle
=
|\mathrm{right}\rangle
$$

### Left

$$
|\psi_c\rangle
=
|\mathrm{left}\rangle
$$

### Symmetric

$$
|\psi_c\rangle
=
\frac{1}{\sqrt{2}}
\left(
|\mathrm{right}\rangle
+
i|\mathrm{left}\rangle
\right).
$$

The symmetric state is useful for producing a spatially symmetric Hadamard walk.

---

## Decoherence

Decoherence is introduced through **Kraus operators** acting on the coin space.

For a channel with operators \(K_k\), the density-operator evolution is

$$
\rho
\longmapsto
\sum_k K_k \rho K_k^\dagger.
$$

This project uses a trajectory-based implementation: at each step a Kraus outcome is sampled according to its quantum probability, the corresponding operator is applied, and the resulting state is normalized.

---

## Dephasing Channel

The dephasing channel suppresses phase coherence while preserving the preferred coin-basis populations.

The implemented Kraus operators are

$$
K_0 = \sqrt{1-p}\,I,
$$

$$
K_1 =
\sqrt{p}\,
|\mathrm{right}\rangle
\langle\mathrm{right}|,
$$

$$
K_2 =
\sqrt{p}\,
|\mathrm{left}\rangle
\langle\mathrm{left}|.
$$

Here \(p\in[0,1]\) is the decoherence strength.

As \(p\) increases, interference is progressively suppressed.

---

## Depolarizing Channel

The depolarizing channel introduces stronger mixing of the coin state using the Pauli operators.

The Kraus operators are

$$
K_0 = \sqrt{1-p}\,I,
$$

$$
K_1 = \sqrt{\frac{p}{3}}X,
\qquad
K_2 = \sqrt{\frac{p}{3}}Y,
\qquad
K_3 = \sqrt{\frac{p}{3}}Z.
$$

The Pauli matrices are

$$
X=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix},
\qquad
Y=
\begin{pmatrix}
0&-i\\
i&0
\end{pmatrix},
\qquad
Z=
\begin{pmatrix}
1&0\\
0&-1
\end{pmatrix}.
$$

---

## Monte Carlo Quantum Trajectories

Instead of evolving a full density matrix, the noisy quantum walk is simulated using independent stochastic quantum trajectories.

For trajectory \(j\), the position distribution at time \(t\) is

$$
P_t^{(j)}(x)
=
|a_t^{(j)}(x)|^2
+
|b_t^{(j)}(x)|^2.
$$

The Monte Carlo estimate of the noisy distribution is

$$
\widehat{P}_t(x)
=
\frac{1}{N}
\sum_{j=1}^{N}
P_t^{(j)}(x),
$$

where \(N\) is the number of simulated trajectories.

Increasing \(N\) generally reduces Monte Carlo variability at the cost of additional computation.

---

## Classical Baseline

The project also computes the exact distribution of a symmetric classical random walk.

At each step,

$$
P(\text{left}) = P(\text{right}) = \frac{1}{2}.
$$

After \(t\) steps, the probability of being at position \(x\) is

$$
P(X_t=x)
=
\binom{t}{(t+x)/2}
\left(\frac12\right)^t,
$$

whenever \(t+x\) is even and \(|x|\le t\).

For this model,

$$
E[X_t]=0,
\qquad
\operatorname{Var}(X_t)=t.
$$

This gives a natural benchmark for identifying the quantum-to-classical transition caused by decoherence.

---

## What the Project Compares

For each simulation, the project can compare:

| Quantity | Noisy QRW | Noiseless QRW | Classical RW |
|---|---:|---:|---:|
| Final position distribution | Yes | Yes | Yes |
| Mean position | Yes | Yes | Yes |
| Variance | Yes | Yes | Yes |
| Variance over time | Yes | Yes | Yes |
| Decoherence channel | Yes | No | No |
| Monte Carlo trajectories | Yes | No | No |

A central question is whether the noisy quantum walk becomes increasingly classical as the decoherence strength \(p\) increases.

---

## R Package

The repository is structured as an R package named `decoherentqrw`.

### Requirements

- R >= 3.5
- `tidyverse`
- `mvtnorm`
- `plotly`

For development and testing:

- `testthat`
- `knitr`
- `rmarkdown`

---

## Installation

Install the development version directly from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "Mahan-H/Decoherence-in-a-Hadamard-Quantum-Random-Walk-on-a-Line"
)
```

Then load the package:

```r
library(decoherentqrw)
```

Alternatively, clone the repository:

```bash
git clone https://github.com/Mahan-H/Decoherence-in-a-Hadamard-Quantum-Random-Walk-on-a-Line.git
cd Decoherence-in-a-Hadamard-Quantum-Random-Walk-on-a-Line
```

---

## Quick Start

### Noiseless Hadamard quantum walk

```r
library(decoherentqrw)

noiseless <- sim_noiseless_qrw(
  T = 20,
  init_coin = "symmetric"
)

noiseless$means
noiseless$vars
```

---

### Noisy quantum walk with dephasing

```r
noisy <- sim_noisy_qrw(
  T = 20,
  N = 500,
  channel = "dephasing",
  init_coin = "symmetric",
  p = 0.10,
  seed = 1
)

noisy$means
noisy$vars
```

---

### Noisy quantum walk with depolarizing noise

```r
noisy_depolarized <- sim_noisy_qrw(
  T = 20,
  N = 500,
  channel = "depolarizing",
  init_coin = "symmetric",
  p = 0.10,
  seed = 1
)
```

---

### Classical symmetric random walk

```r
classical <- sim_srw(T = 20)

classical$means
classical$vars
```

---

## Compare All Three Models

```r
T <- 20

noisy <- sim_noisy_qrw(
  T = T,
  N = 500,
  channel = "dephasing",
  init_coin = "symmetric",
  p = 0.10,
  seed = 1
)

noiseless <- sim_noiseless_qrw(
  T = T,
  init_coin = "symmetric"
)

classical <- sim_srw(T)

final_distribution <- build_dist(
  noisy,
  noiseless,
  classical
)

variance_comparison <- build_variance_overlay(
  noisy,
  noiseless,
  classical
)

summary_table <- build_summary_table(
  noisy,
  noiseless,
  classical
)

head(final_distribution)
head(variance_comparison)
summary_table
```

---

## Main Functions

### Quantum objects and state construction

| Function | Purpose |
|---|---|
| `pos_grid()` | Construct the position grid |
| `hadamard_coin()` | Construct the Hadamard coin operator |
| `pauli_X()` | Construct the Pauli \(X\) operator |
| `pauli_Y()` | Construct the Pauli \(Y\) operator |
| `pauli_Z()` | Construct the Pauli \(Z\) operator |
| `basis_right()` | Right coin basis vector |
| `basis_left()` | Left coin basis vector |
| `proj_right()` | Projector onto the right state |
| `proj_left()` | Projector onto the left state |
| `identity_coin()` | Coin-space identity |
| `init_coin_state()` | Construct an initial coin state |
| `initialize_state()` | Initialize the full quantum-walk state |

### Evolution and decoherence

| Function | Purpose |
|---|---|
| `apply_coin_op()` | Apply a coin-space operator |
| `apply_shift_op()` | Apply the conditional position shift |
| `apply_unitary_step()` | Perform one noiseless Hadamard step |
| `kraus_dephasing()` | Construct dephasing Kraus operators |
| `kraus_depolarizing()` | Construct depolarizing Kraus operators |
| `get_kraus_ops()` | Select the requested decoherence channel |
| `kraus_probs()` | Compute Kraus-outcome probabilities |
| `apply_decoherence_step()` | Sample and apply one decoherence event |

### Simulations and summaries

| Function | Purpose |
|---|---|
| `sim_noisy_qrw()` | Simulate a noisy QRW using Monte Carlo trajectories |
| `sim_noiseless_qrw()` | Simulate the coherent Hadamard QRW |
| `sim_srw()` | Compute the classical symmetric random walk |
| `pos_dist()` | Compute a position probability distribution |
| `mean_pos()` | Compute mean position |
| `var_pos()` | Compute position variance |
| `build_dist()` | Compare final distributions |
| `build_variance_overlay()` | Compare variance over time |
| `build_summary_table()` | Build a final summary table |

---

## Simulation Parameters

### `T`

Number of quantum-walk steps.

A walk of \(T\) steps is represented on the position grid

$$
-T,-T+1,\ldots,T-1,T.
$$

### `N`

Number of Monte Carlo quantum trajectories used for the noisy walk.

A larger value generally provides a more stable approximation to the noisy probability distribution.

### `channel`

Supported values:

```r
"dephasing"
"depolarizing"
```

### `init_coin`

Supported values:

```r
"right"
"left"
"symmetric"
```

### `p`

Decoherence strength.

```text
0 <= p <= 1
```

At \(p=0\), the selected channel introduces no noise.

### `seed`

Optional random seed for reproducible Monte Carlo simulations.

---

## Interactive Shiny Application

The repository includes a Shiny interface under `shiny-app/`.

The application allows users to control:

- number of walk steps \(T\)
- number of noisy trajectories \(N\)
- initial coin state
- decoherence channel
- decoherence strength \(p\)
- random seed
- which models are displayed
- overlay or separate plotting modes

The app provides four main views:

### Final Distribution

Compares the probability of the walker occupying each position after the final time step.

Quantum interference typically produces sharper oscillatory structures, while stronger decoherence tends to smooth the distribution.

### Variance vs. Time

Shows how quickly each model spreads.

This is useful for distinguishing:

- ballistic quantum spreading
- intermediate decohering behaviour
- diffusive classical spreading

### Final Summary

Reports final-time mean and variance for the selected models.

### Model Guide

Provides an interactive explanation of the mathematical models and interpretation of the plots.

---

## Running the Shiny App Locally

Install the required packages:

```r
install.packages(
  c(
    "shiny",
    "bslib",
    "ggplot2",
    "tidyr",
    "plotly",
    "tidyverse",
    "mvtnorm"
  )
)
```

Then, from the project root:

```r
shiny::runApp("shiny-app")
```

---

## Interpreting the Results

### Noiseless quantum walk

A coherent Hadamard walk typically exhibits:

- interference
- oscillatory probability peaks
- greater probability mass away from the origin
- ballistic spreading

### Classical random walk

The classical symmetric random walk typically exhibits:

- no quantum interference
- a smoother distribution
- concentration closer to the origin
- variance growing linearly with time

### Noisy quantum walk

The noisy walk lies between these two regimes.

As decoherence becomes stronger, quantum interference is increasingly suppressed. A key empirical signature of the quantum-to-classical transition is that the noisy walk's distribution and variance begin to resemble those of the classical random walk.

---

## Example Research Questions

The code can be used to investigate questions such as:

1. How does dephasing alter the final position distribution?
2. How does depolarizing noise compare with dephasing?
3. At what decoherence strength does the quantum walk begin to resemble a classical walk?
4. How does the variance scaling change as noise increases?
5. How sensitive is the walk to the initial coin state?
6. How many Monte Carlo trajectories are required for stable numerical estimates?

---

## Testing

The project uses `testthat`.

Run the test suite with:

```r
devtools::test()
```

or:

```r
testthat::test_dir("tests/testthat")
```

The tests check important mathematical and numerical properties, including:

- the position grid
- Hadamard unitarity
- Pauli-matrix identities
- basis vectors and projectors
- normalization of quantum states
- conditional shifts
- Kraus completeness
- Kraus probabilities
- noisy-state normalization
- exact small classical distributions
- probability normalization
- expected output structures

---

## Repository Structure

```text
.
├── R/
│   └── main.R
├── man/
├── shiny-app/
│   ├── R/
│   ├── app.R
│   ├── server-plots.R
│   └── manifest.json
├── tests/
│   ├── testthat/
│   │   └── tests.R
│   └── testthat.R
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── LICENSE.md
├── STA380-Project.Rproj
└── README.md
```

---

## Reproducibility

The noisy simulation is stochastic because Kraus outcomes are sampled during each quantum trajectory.

For reproducible results, specify a seed:

```r
result <- sim_noisy_qrw(
  T = 20,
  N = 500,
  channel = "dephasing",
  init_coin = "symmetric",
  p = 0.1,
  seed = 123
)
```

Using the same simulation settings and seed should reproduce the same Monte Carlo result.

---

## Computational Considerations

The computational cost of the noisy model increases with both:

- the number of steps \(T\)
- the number of trajectories \(N\)

Larger values of `N` reduce Monte Carlo variability but increase runtime.

For interactive exploration, moderate values of `T` and `N` are usually sufficient. Larger values are more appropriate when estimating distributions or variance curves for formal analysis.

---

## Motivation

Quantum random walks are quantum analogues of classical random walks and are important in quantum information, quantum algorithms, transport models, and the study of quantum-to-classical transitions.

The central idea of this project is that the difference between quantum and classical spreading is not simply a different probability rule. It results from **coherent interference between quantum amplitudes**.

Decoherence progressively destroys that interference.

By simulating the coherent walk, noisy walk, and classical walk in the same framework, this project provides a computational demonstration of how environmental noise can transform quantum dynamics into increasingly classical behaviour.

---

## Author

**Mahan Hooshmandkhayat**

University of Toronto

---

## License

This project is distributed under the **MIT License**. See `LICENSE` and `LICENSE.md` for details.
