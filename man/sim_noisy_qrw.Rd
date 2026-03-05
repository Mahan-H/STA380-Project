\name{sim_noisy_qrw}
\alias{sim_noisy_qrw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_noisy_qrw(T, N, channel, init_coin, p, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{channel}{
%%     ~~Describe \code{channel} here~~
}
  \item{init_coin}{
%%     ~~Describe \code{init_coin} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (T, N, channel, init_coin, p, seed = NULL) 
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    pos <- pos_grid(T)
    times <- 0:T
    dists <- matrix(0, nrow = length(times), ncol = length(pos))
    for (i in 1:N) {
        state <- initialize_state(T, init_coin)
        dists[1, ] <- dists[1, ] + pos_dist(state)
        for (t in 1:T) {
            state <- apply_unitary_step(state)
            state <- apply_decoherence_step(state, channel, p)
            dists[t + 1, ] <- dists[t + 1, ] + pos_dist(state)
        }
    }
    normalize_dists <- dists/N
    means <- numeric(length(times))
    vars <- numeric(length(times))
    for (t in 0:T) {
        probs <- normalize_dists[t + 1, ]
        means[t + 1] <- mean_pos(pos, probs)
        vars[t + 1] <- var_pos(pos, probs)
    }
    data <- list(T = T, dists = normalize_dists, means = means, 
        vars = vars)
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
