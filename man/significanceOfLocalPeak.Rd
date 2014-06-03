\name{significanceOfLocalPeak}
\alias{significanceOfLocalPeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical significance of local periodogram peaks
}
\description{
Calculate statistical significance for a secondary periodogram peak (i.e. a non-global periodogram maximum), based on the null hypothesis of an OUSS process.
}
\usage{
significanceOfLocalPeak(power_o, lambda, power_e, time_step, 
                        Nfreq, peakFreq, peakPower)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power_o}{
    Positive number. Power at zero-frequency stemming from the underlying OU process.
}
  \item{lambda}{
    Positive number. Resilience (or relaxation rate) of the OU process, in inverse time units. This is also the inverse correlation time of the OU process. 
}
  \item{power_e}{
    Non-negative number. Asymptotic power at large frequencies due to random measurement errors. Setting this to zero corresponds to the classical OU process.
}
  \item{time_step}{
    Positive number. The time step of the time series that was used to calculate the periodogram.
}
  \item{Nfreq}{
    The number of frequencies from which the local periodogram peak was picked. Typically equal to the number of frequencies in the periodogram.
}
  \item{peakFreq}{
    Single number. The frequency of the focal peak.
}
  \item{peakPower}{
    Single number. The periodogram power calculated for the focal peak.
}
}
\details{
The OUSS parameters \code{power_o}, \code{lambda} and \code{power_e} will typically be maximum-likelihood fitted values returned by \code{\link{evaluate.pm}}. \code{time_step} is also returned by \code{evaluate.pm} and is inferred from the analysed time series. The examined periodogram peak (as defined by \code{peakFreq}) will typically be a secondary peak of interest, masked by other stronger peaks or a low-frequency maximum. The significance of such a peak is not defined by standard tests.
}
\value{
The returned P-value (referred to as ``local P-value'') is the probability that an OUSS process with the specified parameters would generate a periodogram with a power-to-expectation ratio greater than \code{peakPower/E}, where \code{E} is the power spectrum of the OUSS process at frequency \code{peakFreq}. Hence, the significance is a measure for how much the peak power deviates from its expectation. The calculated value is an approximation. It becomes exact for long regular time series.
}
\references{
Louca, S., Doebeli, M. - Detecting cyclicity in ecological time series (in review, as of June 2014)
}
\author{
Stilianos Louca
}
\note{
This statistical significance is not equivalent to the one calculated by \code{\link{evaluate.pm}} for the global periodogram maximum.
If the investigated periodogram peak is a global maximum, then the P-value returned by \code{\link{evaluate.pm}} should be preferred, as it also takes into account the absolute magnitude of the peak.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{evaluate.pm}}
}
\examples{
  # In this example we generate a random cyclic time series, where the peak is (most likely)
  # masked by a strong low-frequency maximum.
  # We will use significanceOfLocalPeak() to evaluate its significance
  # based on its deviation from the expected power.
	
  # generate cyclic time series by adding a periodic signal to an OUSS process
  period    = 1;
  times     = seq(0,20,0.2);
  signal    = 0.75 * cos(2*pi*times/period) +
              generate_ouss(times, mu=0, sigma=1, lambda=1, epsilon=0.5);
  
  # calculate periodogram and fit OUSS model
  report    = evaluate.pm(times=times, signal=signal);
  print(report)
	
  # find which periodogram mode approximately corresponds to the frequency we are interested in
  cycleMode = which(report$frequencies>=0.99/period)[1]; 
	
  # calculate P-value for local peak
  Pvalue    = significanceOfLocalPeak(power_o   = report$power_o, 
                                      lambda    = report$lambda, 
                                      power_e   = report$power_e, 
                                      time_step = report$time_step, 
                                      Nfreq     = length(report$frequencies), 
                                      peakFreq  = report$frequencies[cycleMode], 
                                      peakPower = report$periodogram[cycleMode]);

  # plot time series
  old.par <- par(mfrow=c(1, 2));
  plot(ts(times), ts(signal), 
       xy.label=FALSE, type="l", 
       ylab="signal", xlab="time", main="Time series (cyclic)", 
       cex=0.8, cex.main=0.9);

  # plot periodogram
  title = sprintf("Periodogram OUSS analysis\nfocusing on local peak at freq=\%.3g\nPlocal=\%.2g",
                  report$frequencies[cycleMode],Pvalue);
  plot(ts(report$frequencies), ts(report$periodogram), 
       xy.label=FALSE, type="l", 
       ylab="power", xlab="frequency", main=title, 
       col="black", cex=0.8, cex.main=0.9);
	
  # plot fitted OUSS power spectrum
  lines(report$frequencies, report$fittedPS, col="red");
	par(old.par)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
