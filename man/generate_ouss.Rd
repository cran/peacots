\name{generate_ouss}
\alias{generate_ouss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate random time series of the OUSS process
}
\description{
Generate a random time series of the 1-dimensional stationary Ornstein-Uhlenbeck state space (OUSS) process.
}
\usage{
generate_ouss(times, mu, sigma, lambda, epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
	Numeric vector of times for which to evaluate OUSS model. Times need to be strictly ascending.
}
  \item{mu}{
	Single number. Deterministic equilibrium of OU process, i.e., the expected value of the time series at any particular time.
}
  \item{sigma}{
	Single number. Standard deviation of OU fluctuations around equilibrium.
}
  \item{lambda}{
	Single non-negative number. Resilience (also known as relaxation rate) of the OU process. This is the inverse of the OU correlation time.
}
  \item{epsilon}{
	Single number. Standard deviation of Gaussian measurement error. Setting this to zero will yield a time series from the classical OU process.
}
}
\details{
The OUSS model describes the measurement of an Ornstein-Uhlenbeck (OU) stochastic process at discrete times with additional uncorrelated Gaussian measurement errors. The OU process itself is a continuous-time random walk (Brownian motion) with linear stabilizing forces, described by the stochastic differential equation
\deqn{
	dX = \lambda(\mu-X) dt + s dW,
}
where \eqn{W} is the standard Wiener process and \eqn{s^2=2\lambda\sigma^2.} The OUSS model is obtained by adding uncorrelated Gaussian numbers with zero mean and variance \eqn{\epsilon^2} to the time series.
}
\value{
A numeric vector of same length as \code{times}, containing sampled values of the OUSS process. These values will all have the same expectation (\code{mu}) and variance (\code{sigma^2+epsilon^2}) but will be correlated.
}
\references{
Louca, S., Doebeli, M. - Detecting cyclicity in ecological time series (in review, as of June 2014)

Dennis, B., Ponciano, J.M. - Density dependent state-space model for population abundance data with unequal time intervals, Ecology (in press as of June 2014)
}
\author{
Stilianos Louca
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ps_ouss}}
}
\examples{
# define times
times = seq(0,100,0.5);

# generate OUSS time series
signal 	= generate_ouss(times=times, mu=0, sigma=1, lambda=1, epsilon=0.5);

# plot time series
plot(ts(times), ts(signal), 
     xy.label=FALSE, type="l", 
     ylab="signal", xlab="time", main="OUSS time series");
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics }
\keyword{ Mathematics }% __ONLY ONE__ keyword per line
