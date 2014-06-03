\name{peacots-package}
\alias{peacots-package}
\alias{peacots}
\docType{package}
\title{
Periodogram peaks in correlated time series
}
\description{
Calculate the periodogram of a time series, maximum-likelihood fit an Ornstein-Uhlenbeck state space (OUSS) null model to the periodogram and evaluate the statistical significance of periodogram peaks against the OUSS null hypothesis.
}
\details{
\tabular{ll}{
Package: \tab peacots\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-03\cr
License: \tab GPL-3\cr
}
The OUSS is a parsimonious model for a stochastically fluctuating variable (e.g. population size) with linear stabilizing forces, subject to uncorrelated measurement errors. Periodogram peaks (putative periodicities) are evaluated against this null hypothesis. In contrast to the white noise null model (the classical null model against which cyclicity is often evaluated), the OUSS process accounts for non-zero correlations between measurements and corrects for the resulting increased power at low frequencies.

Use \code{\link{evaluate.pm}} to calculate the periodogram of a time series, fit the OUSS null model and calculate the statistical significance of the periodogram maximum.

Use \code{\link{plotReport}} to generate a simple plot of the results returned by \code{\link{evaluate.pm}}.

Use \code{\link{significanceOfLocalPeak}} to evaluate the statistical significance of a secondary peak (i.e. non-global maximum) in the periodogram.

Use \code{\link{runExample}} to run an example \code{peacots} analysis based on simulation data.

Use \code{\link{evaluate.pm.wn}} to evaluate the statistical significance of the periodogram maximum against the white noise null hypothesis. This is the classical test, included for comparison.

Use \code{\link{ps_ouss}} to calculate the power spectrum of a particular OUSS process.

Use \code{\link{generate_ouss}} to generate random time series of a particular OUSS process.
}
\author{
Stilianos Louca
}
\references{
Louca, S., Doebeli, M. - Detecting cyclicity in ecological time series (in review, as of June 2014)
}
\keyword{ package }
\keyword{ statistics }
\examples{
# Generate a cyclic time series and analyse using peacots

# Parameters
lambda         = 1;    # inverse correlation time of OU process
cyclePeriod    = 1;
cycleAmplitude = 0.75;
times          = seq(0,20,0.2);

# Example 1

# generate cyclic time series by adding a periodic signal to an OUSS process
signal 	= cycleAmplitude * cos(2*pi*times/cyclePeriod) +
          generate_ouss(times, mu=0, sigma=1, lambda=lambda, epsilon=0.5);

# Find periodogram peak and estimate statistical significance
# Ignore frequencies lower than a pre-defined threshold
# to avoid masking by low-frequency maximum
report 	= evaluate.pm(times=times, signal=signal, 
                      minPeakFreq=lambda/3, 
                      minFitFreq=lambda/3);

# plot overview of periodogram peak analysis
plotReport(sprintf("Cyclic at frequency \%.2g",1/cyclePeriod), 
           times=times, signal=signal, report=report);


# Example 2 (using the same time series)
# In this example we don't use low-frequency trimming
# Instead, we will focus on a particular (local) periodogram peak
# and estimate its 'local' statistical significance

# calculate periodogram and fit OUSS model
report    = evaluate.pm(times=times, signal=signal);

# find the periodogram mode approximately corresponding to the frequency we are interested in
cycleMode = which(report$frequencies>=0.99/cyclePeriod)[1];

# calculate local P-value for this peak
Pvalue    = significanceOfLocalPeak(power_o   = report$power_o, 
                                    lambda    = report$lambda, 
                                    power_e   = report$power_e, 
                                    time_step = report$time_step, 
                                    Nfreq     = length(report$frequencies), 
                                    peakFreq  = report$frequencies[cycleMode], 
                                    peakPower = report$periodogram[cycleMode]);

# print result
cat(sprintf("Local P-value = \%.3g for peak at frequency=\%.3g\n", 
    Pvalue, report$frequencies[cycleMode]));
}
