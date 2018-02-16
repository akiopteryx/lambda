\name{sim.coord}
\alias{sim.coord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate coordinate data.
}
\description{
Generates 2-D and 3-D data based on a simple variance-covariance matrix.
}
\usage{
sim.coord(n.obj, n.landmark, var=1, cov=0.5, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.obj}{
Number of specimens.
}
  \item{n.landmark}{
Number of landmarks.
}
  \item{var} {
  Single variance applied to each landmark (default=1.0).
  }
  \item{cov} {
  Single covariance applied to each pair of landmarks (default=0.5).
  }
  \item{dim} {
  Number of spatial dimensions (e.g., 2 for 2-D data, 3 for 3-D data)
  }
}
\details{
}
\value{
The function returns a matrix of simulated coordinate values where rows and columns represent specimens and shape variables, respectively.
}
\references{
Watanabe, A. How many landmarks are enough to characterize shape and size variation? ### UPDATE WHEN PUBLISHED ##
}
\author{
Akinobu Watanabe
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sim.3d.data <- sim.coord(n.obj=10, n.landmark=20, var=1, cov=0.5, dim=3)
sim.3d.data
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data}% use one of  RShowDoc("KEYWORDS")
\keyword{simulation}% __ONLY ONE__ keyword per line
