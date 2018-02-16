\name{lasec}
\alias{lasec}
\title{
Landmark Sampling Evaluation Curve
}
\description{
Generate a sampling curve to assess whether sufficient number of landmarks has been sampled.
}
\usage{
lasec(coord.data, n.dim, iter = 1000, show.progress = T)
}
\arguments{
  \item{coord.data}{
2-D or 3-D coordinate data arranged in n x p matrix with n=specimens, p=shape variables. NB: data should be unaligned for centroid size analysis).
}
  \item{n.dim}{
Number of spatial dimensions (i.e., 2 or 3 for 2-D or 3-D).
}
  \item{iter}{
Number of subsampling rounds.
}
}
\details{
The function reads in landmark data without any missing data and generates two sampling curves. "LaSEC_PSSCurve_Shape.pdf" displays how sampling additional landmarks converges to the distribution of shape variation in the parent dataset. "LaSEC_PSSCurve_CS.pdf" displays how sampling additional landmarks converges to the distribution of centroid size variation in the parent dataset. These sampling curves allows the user to determine whether landmark sampling is adequate to capture the shape and size variation of specimens by looking at the extent of the plateau in the curve and constricting fit (1-PSS) values.
}
\value{
Function returns the following:
\item{pss}{A matrix of 1–PSS values where each row comprises sequence of fit values at each subsampling of landmarks with respect to shape variation
}
\item{median.pss}{A vector comprising median 1–PSS values at each subsampling of landmarks
}
\item{pss.cs}{A matrix of 1–PSS values where each row comprises sequence of fit values at each subsampling of landmarks with respect to centroid size variation
}
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
data(crocskull)
lasec(coord.data=crocskull, n.dim=3, iter=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sampling}% use one of  RShowDoc("KEYWORDS")
\keyword{data assessment}% __ONLY ONE__ keyword per line
