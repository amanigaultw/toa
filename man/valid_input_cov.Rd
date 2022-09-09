\name{valid_input_cov}
\alias{valid_input_cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
check if cov matrix/vector is valid
}
\description{
valid_input_cov checks whether the cov matrix/vector is valid for use by regressions that model each gene individually, including checking for missing values and checking if any covariate variables are a constants.
}
\usage{
valid_input_cov(cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov}{
a matrix or vector of variables to control for when examining differential gene expression as a function of the x variable.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
}
