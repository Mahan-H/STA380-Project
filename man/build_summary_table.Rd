\name{build_summary_table}
\alias{build_summary_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
build_summary_table(noisy_w, noiseless_w, srw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noisy_w}{
%%     ~~Describe \code{noisy_w} here~~
}
  \item{noiseless_w}{
%%     ~~Describe \code{noiseless_w} here~~
}
  \item{srw}{
%%     ~~Describe \code{srw} here~~
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
function (noisy_w, noiseless_w, srw) 
{
    T <- noisy_w$T
    noisy_mean <- noisy_w$means[T + 1]
    noisy_var <- noisy_w$vars[T + 1]
    noiseless_mean <- noiseless_w$means[T + 1]
    noiseless_var <- noiseless_w$vars[T + 1]
    srw_mean <- srw$means[T + 1]
    srw_var <- srw$vars[T + 1]
    summary <- data.frame(c("noisy", "noiseless", "classical"), 
        c(noisy_mean, noiseless_mean, srw_mean), c(noisy_var, 
            noiseless_var, srw_var))
    names(summary)[1] <- "model"
    names(summary)[2] <- "mean"
    names(summary)[3] <- "variance"
    return(summary)
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
