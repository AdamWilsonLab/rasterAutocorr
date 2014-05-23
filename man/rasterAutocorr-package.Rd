\name{rasterAutocorr-package}
\alias{rasterAutocorr-package}
\alias{rasterAutocorr}
\docType{package}
\title{
Quickly calculate spatial autocorrelation on 2D rasters
}
\description{
Calculates spatial autocorrelation on a two-dimensional raster* object using fft() methods.
Install with: \cr
\code{ \cr
library(devtools)  \cr
install_github("adammwilson/rasterAutocorr")}
}
\details{
\tabular{ll}{
Package: \tab rasterAutocorr\cr
Type: \tab Package\cr
Version: \tab 0.0001\cr
Date: \tab 2014-05-23\cr
License: \tab \url{http://www.r-project.org/Licenses/GPL-3}\cr
}
\code{autocor} calculates spatial autocorrelation.
}
\author{
Adam M. Wilson

Maintainer: Who to complain to Adam M. Wilson <adam.wilson@yale.edu>
}
\references{
Xianlin Ma, Tingting Yao, A program for 2D modeling (cross) correlogram tables using fast Fourier transform, Computers & Geosciences, Volume 27, Issue 7, August 2001, Pages 763-774, ISSN 0098-3004, \url{http://dx.doi.org/10.1016/S0098-3004(01)00007-3}.

Autocorrelation and Self Convolution \url{http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm}

\url{http://octave.sourceforge.net/signal/function/xcorr2.html}

\url{http://www.mathworks.com/help/matlab/ref/fftshift.html}
}
\keyword{ package }
\seealso{
\code{\link[raster:raster-package]{raster}}
\code{\link[fft:base]{fft}}
}
%\examples{
%}
