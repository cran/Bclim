\name{NIGB}
\alias{NIGB}
\title{
Normal-Inverse Gaussian bridging function
}
\description{
Normal-Inverse Gaussian bridging function
}
\usage{
NIGB(delta, IG.sum, tb.points, c.start, c.end)
}
\arguments{
  \item{delta}{
IG parameter 1
}
  \item{IG.sum}{
Sum of IG values
}
  \item{tb.points}{
Time points over which to interpolate
}
  \item{c.start}{
Climate start value
}
  \item{c.end}{
Climate end value
}
}
\details{
Not meant for use outside of \code{\link{BclimInterp}}.
}
\value{
Interpolated NIG values returned to \code{\link{BclimInterp}}.
}
\author{
Thinh K Doan
}
\seealso{
\code{\link{Bclim}}
}
\examples{
## See other Bclim functions
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
