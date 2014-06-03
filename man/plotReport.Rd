\name{plotReport}
\alias{plotReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results of a \code{evaluate.pm} analysis
}
\description{
Create a simple plot of a time series and the results of a \code{\link{evaluate.pm}} analysis (including the periodogram and the fitted OUSS power spectrum).
}
\usage{
plotReport(name="", times=NULL, signal=NULL, 
           report=NULL, plotFile=NULL, dataFile=NULL, sep=" ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
    Character. A short name for the time series to be used for the plots (e.g. 'long-term study' or 'hare population').
}
  \item{times}{
    Numeric vector. The time points of the time series used for the analysis. Set to \code{NULL} to not plot the time series.
}
  \item{signal}{
    Numeric vector. The time series values (signal) used for the analysis. Set to \code{NULL} to not plot the time series.
}
  \item{report}{
    The value returned by \code{evaluate.pm}. This is a list with the periodogram analysis results. Set to \code{NULL} to not plot any analysis results.
}
  \item{plotFile}{
    An optional path to a PDF file to be generated with the plot.
}
  \item{dataFile}{
    An optional path to a data file for storing the time series and the results of the analysis.
}
  \item{sep}{
    Separator to be used for the data file. Only relevant if \code{dataFile} is not \code{NULL}.
}
}
\references{
Louca, S., Doebeli, M. - Detecting cyclicity in ecological time series (in review, as of June 2014)
}
\author{
Stilianos Louca
}
\seealso{
\code{\link{evaluate.pm}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# generate cyclic time series by adding a periodic signal to an OUSS process
times  = seq(0,20,0.2);
signal = 0.75 * cos(2*pi*times) + generate_ouss(times, mu=0, sigma=1, lambda=1, epsilon=0.5);

# find periodogram peak and estimate statistical significance
report = evaluate.pm(times=times, signal=signal, minPeakFreq=0.3, minFitFreq=0.3);

# plot overview of periodogram peak analysis
plotReport(sprintf("Example"), times=times, signal=signal, report=report);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ statistics }
