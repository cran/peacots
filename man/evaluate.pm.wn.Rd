\name{evaluate.pm.wn}
\alias{evaluate.pm.wn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical significance of periodogram peaks (classical)
}
\description{
Calculates the Lomb-Scargle periodogram for the given time series and estimates the statistical significance of the global periodogram maximum based on the null hypothesis of uncorrelated (white) noise. Available for historical reasons and for comparison purposes.
}
\usage{
evaluate.pm.wn(times, signal, minPeakFreq=0, minFitFreq=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
     Numerical vector. Time points of the time series.
}
  \item{signal}{
     Numerical vector of same size as \code{times}. Values of the time series.
}
  \item{minPeakFreq}{
    Single non-negative number. Minimum considered frequency when determining periodogram peak. Use this to ignore low-frequency components from the spectrum.
}
  \item{minFitFreq}{
    Single non-negative number. Minimum considered frequency when fitting the white noise null model to the periodogram. Use this to ignore low-frequency components from the spectrum.
}
}
\value{
A list with the entries
\item{error}{Will be \code{TRUE} if an error occured, \code{FALSE} otherwise.}
\item{errorMessage}{A short error message if \code{error==TRUE}.}
If \code{error==FALSE}, the returned list also includes:
\item{frequencies}{Available periodogram frequencies as a numerical vector.}
\item{periodogram}{Periodogram powers corresponding to the returned frequencies, as a numerical vector.}
\item{peakMode}{An integer indicating the position of the global periodogram maximum (starting at \code{minPeakFreq}) in the vector \code{frequencies}.}
\item{powerEstimate}{The estimated white noise power. Estimated from the average periodogram power, which corresponds to using the total variance of the time series (if \code{minFitFreq==0}).}
\item{minPeakMode}{The minimum periodogram mode considered for determining the periodogram peak. This will be \code{1} if \code{minPeakFreq==0}.}
\item{minFitMode}{The minimum periodogram mode considered for estimating the white noise power. This will be \code{1} if \code{minFitFreq==0}.}
\item{RSS}{The sum of squared residuals of the periodogram from the estimated white noise power.}
\item{P}{Statistical significance of periodogram peak. This is the probability that a white noise periodogram (of the estimated power) would generate a peak at least as strong as the observed peak (among the considered frequencies). The calculated P-value is only an approximation that becomes exact for long regular time series.}
}
\references{
Scargle, J. D. (1982) - Studies in astronomical time series analysis. II Statistical aspects of spectral analysis of unevenly spaced data, The Astrophysical Journal 263, 835--853

Horne, J. H., Baliunas, S. L. (1986) - A prescription for period analysis of unevenly sampled time series, The Astrophysical Journal 302, 757--763

Louca, S., Doebeli, M. (2015) Detecting cyclicity in ecological time series, Ecology 96: 1724--1732
}
\author{
Stilianos Louca
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{evaluate.pm}}
}
\examples{
# generate time series
times  = seq(0,20,0.25);
signal = rnorm(n=length(times));
report = evaluate.pm.wn(times=times, signal=signal);

# plot time series
old.par <- par(mfrow=c(1, 2));
plot(ts(times), ts(signal), 
     xy.label=FALSE, type="l", 
     ylab="signal", xlab="time", main="OUSS time series");
     
# plot periodogram
title = sprintf("Periodogram OUSS analysis\n(peak freq=\%.3g, P=\%.2g)",
                report$frequencies[report$peakMode],report$P);
plot(ts(report$frequencies),
     ts(report$periodogram), 
     xy.label=FALSE, type="l", 
     ylab="power", xlab="frequency", 
     main=title, col="black");
	
# plot fitted flat WN power
lines(c(report$frequencies[1],tail(report$frequencies,1)), 
      c(report$powerEstimate, report$powerEstimate ), col="blue");
par(old.par)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ statistics }
