\name{sim_noiseless_qrw}
\alias{sim_noiseless_qrw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_noiseless_qrw(T, init_coin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{init_coin}{
%%     ~~Describe \code{init_coin} here~~
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
function (T, init_coin) 
{
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
    for (t in 1:T) {
        state <- apply_unitary_step(state)
        prob <- pos_dist(state)
        dists[t + 1, ] <- prob
        means[t + 1] <- mean_pos(pos, prob)
        vars[t + 1] <- var_pos(pos, prob)
    }
    data <- list(T = T, dists = dists, means = means, vars = vars)
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
