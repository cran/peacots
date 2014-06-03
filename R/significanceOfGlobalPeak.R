significanceOfGlobalPeak <-
function(power_o, lambda, power_e, time_step, frequencies, peakFreq, peakPower, accuracy=0.005){
	
	#Number of trials for estimation of false alarm probability P.
	#The variance of the estimator will be P(1-P)/trials and at most 0.25/trials
	trials 		= max(10,ceiling(0.25/accuracy^2));
	
	#Preliminary preparations
	countPositives	= 0;
	EpeakPower 		= ps_ouss(peakFreq, power_o, lambda, power_e, time_step);
	Epowers 		= ps_ouss(frequencies, power_o, lambda, power_e, time_step);
	
	#Generate all exponentially distributed variables corresponding to power estimates in a hypothetical periodogram
	rexps = matrix(0, trials, length(frequencies));
	for(m in 1:length(frequencies)){
		if(Epowers[m]==0){
			rexps[,m] = rep(0, trials);
		}else{
			rexps[,m] 	= rexp(n=trials, rate=1/Epowers[m]);
			if(is.na(rexps[1,m])){
				# something went wrong
				cat(sprintf(" REXPS generated NA: Epowers[%a]=%g, power_o=%g, lambda=%g, power_e=%g, time_step=%g\n", m, Epowers[m], power_o, lambda, power_e, time_step));
			}
		}
	}
	
	#Evaluate trial periodograms, keeping track of how many ended up having peaks at least as extreme as the case given
	for(n in 1:trials){
		m = which.max(rexps[n,]);		# detect peak in current trial
		r_peakFreq 	= frequencies[m];
		r_peakPower = rexps[n,m];
		if((r_peakPower>=peakPower) && (r_peakPower/Epowers[m] >= peakPower/EpeakPower)){
			countPositives = countPositives + 1;
		}
	}
	
	#Estimated false alarm probability is fraction of random peridograms with stronger peaks than case given 
	return(countPositives/trials);
}
