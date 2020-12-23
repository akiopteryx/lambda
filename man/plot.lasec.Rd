\name{plot.lasec}
\alias{plot.lasec}
\title{
Plot results of LaSEC analysis.
}
\description{
Generates landmark sampling curves from LaSEC results.
}
\usage{
plot.lasec(lasec.result)
}
\arguments{
  \item{lasec.result}{
Variable recording entire LaSEC output.
}
}
\details{
The function allows user to generate landmark sampling curves from LaSEC results on their own.
}
\value{
The function creates a landmark sampling curve to assess the quality of shape data.
}
\references{
Watanabe, A. 2018. How many landmarks are enough to characterize shape and size variation? PLoS ONE 13(6): e0198341. https://doi.org/10.1371/journal.pone.0198341
}
\author{
Akinobu Watanabe
}

\examples{
data(crocskull)
lasec.results <- lasec(coord.data=crocskull, n.dim=3, iter=10)
plot.lasec(lasec.results)
}

\keyword{data}
\keyword{simulation}
