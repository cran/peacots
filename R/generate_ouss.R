generate_ouss <-
function(	times,		# time points
			mu,			# deterministic equilibrium
			sigma,		# standard deviation of OU process (i.e. of its stationary distribution)
			lambda,		# relaxation rate of OU process
			epsilon){	# standard deviation of Gaussian measurement error. Set this to 0 for a classical OU process.
	N = length(times)
	if(N==0) return(c());
	uncorrelated_signal = rnorm(N,mean=0,sd=sigma);
    signal		= rep(0,N);
    signal[1] 	= uncorrelated_signal[1]; # random starting point
    
    # simulate OU time series through correlated draws
    for(k in 2:N){
		rho	= exp(-lambda*(times[k]-times[k-1])); # correlation between consecutive time steps
        signal[k] = rho*signal[k-1] + sqrt(1-rho^2)*uncorrelated_signal[k];
	}
	
	# shift uniformly
	signal = signal + mu;
	
	# add Gaussian measurement errors to obtain an OU state space model
	signal = signal + rnorm(N,mean=0,sd=epsilon);
	
	return(signal);
}
