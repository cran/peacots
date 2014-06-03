evaluate.pm <-
function(times, signal, minPeakFreq=0, minFitFreq=0, iterations=100, accuracy=0.002, startRadius=2, verbose=FALSE){
	# Pre-calculated FAP correction grid
	# Maps time series length & FAP estimate to CDF
	FAP_correction_grid = 	"20 50 100 200 500,0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1,\
							0 0 0 0 0 0 0 0 0 0 0 0 4.14413e-05 0.000100558 0.000179206 0.000234285 0.000304372 0.000545357 0.000718865 0.00143079 0.00280069 0.00430458 0.00600418 0.00808732 0.0104055 0.0126051 0.0153489 0.0180491 0.0208072 0.0240908 0.0274736 0.0310694 0.0351576 0.0392573 0.0434526 0.0478777 0.0523421 0.103101 0.161997 0.225263 0.293975 0.36512 0.438276 0.508701 0.580871 0.713995 0.839797 0.937171 0.992196 1 \
							2.92994e-05 5.57744e-05 0.000115764 0.000200192 0.000214838 0.000223307 0.000237501 0.00025378 0.000272656 0.000279051 0.000381618 0.000607157 0.000894802 0.00128185 0.00167089 0.00201414 0.00236281 0.00268935 0.00302546 0.00506285 0.00734405 0.0100983 0.0126831 0.0158258 0.0191664 0.0229441 0.0267799 0.0308015 0.0347535 0.0393276 0.044053 0.0489552 0.0538847 0.0595843 0.0646582 0.0700988 0.0759567 0.137677 0.205985 0.275746 0.349259 0.42107 0.493343 0.562025 0.63045 0.75377 0.861418 0.944232 0.990784 1 \
							0.000560907 0.000745208 0.000896268 0.00106066 0.00113213 0.00119029 0.0012824 0.00137423 0.00148657 0.00156486 0.0021311 0.00284119 0.00351032 0.00417341 0.00490792 0.00543324 0.00606696 0.00652527 0.00707935 0.0100837 0.0134424 0.0171353 0.0208599 0.0249107 0.0291754 0.033615 0.0382322 0.0429625 0.0478509 0.0530736 0.0585065 0.0640751 0.0699815 0.0763537 0.0822729 0.0884222 0.0948485 0.162887 0.23569 0.308811 0.383264 0.456772 0.528003 0.595692 0.661056 0.777176 0.875264 0.94882 0.990196 1 \
							0.00276238 0.00316131 0.00342642 0.00366936 0.00381359 0.00397443 0.00416159 0.00430386 0.00445382 0.00459492 0.00561732 0.0066108 0.00765389 0.00853629 0.00942272 0.0102632 0.0110513 0.01169 0.0124388 0.0161743 0.0202517 0.0247817 0.0291928 0.0340112 0.0388258 0.0438408 0.0490462 0.0541426 0.0597469 0.0654283 0.0711345 0.0773415 0.0838424 0.0904408 0.0969427 0.103528 0.110407 0.181655 0.257547 0.332938 0.408315 0.482212 0.553064 0.619748 0.68277 0.794365 0.885726 0.952456 0.990271 1 \
							0.0047867 0.00525947 0.00563846 0.00593054 0.00615252 0.00639303 0.00660491 0.00677807 0.00694337 0.00711613 0.00831689 0.00949392 0.0106534 0.0116221 0.0126563 0.0135687 0.0143976 0.0151618 0.0160323 0.0203308 0.0247287 0.0296046 0.0345213 0.039511 0.0444977 0.0497686 0.0554052 0.0609163 0.0668458 0.0726852 0.0786081 0.0851175 0.0917448 0.0985411 0.105117 0.111844 0.118832 0.192173 0.269398 0.34554 0.421054 0.494716 0.564963 0.630945 0.692703 0.801271 0.888958 0.95224 0.989029 1";

	# average time step
	time_step = (tail(times,n=1)-times[1])/(length(times)-1);
	
	
	#######################
	# Calculate periodogram
	
	
	# Calculate Lomb-Scargle periodogram from time series
	if(verbose) cat(sprintf("Calculating periodogram..\n"));
	spectrum 	= LombScarglePeriodogram(times, signal);
	frequencies	= spectrum$freq;
	powers		= spectrum$spec;
	frequencies = frequencies[-1]; 		#discard zero-frequency mode
	powers	 	= powers[-1]; 			#discard zero-frequency mode
	if(tail(frequencies,n=1)<minPeakFreq) return(list(error=TRUE, errorMessage="All periodogram frequencies are below minPeakFreq"));
	if(tail(frequencies,n=1)<minFitFreq) return(list(error=TRUE, errorMessage="All periodogram frequencies are below minFitFreq"));
	minPeakMode	= which(frequencies>=minPeakFreq)[1];	#minimum mode considered for peak search. Does not influence ML fitting.
	minFitMode	= which(frequencies>=minFitFreq)[1];	#minimum mode considered for maximum-likelihood fitting. Does not influence peak search.
	
	# Detect spectral peak
	peakMode 	= (minPeakMode-1) + which.max(powers[minPeakMode:length(powers)]);
	peakPower 	= powers[peakMode];
	peakFreq 	= frequencies[peakMode];	
	


	######################################
	# Estimate OU process parameters (power_o & lambda & power_e) from periodogram using maximum likelihood
	
	if(verbose) cat(sprintf("Maximum-likelihood estimating OUSS power parameters..\n"));
	mlfit				= mlfit_ouss(	frequencies = frequencies[minFitMode:length(frequencies)], 
										periodogram	= powers[minFitMode:length(powers)], 
										time_step	= time_step, 
										startRadius	= startRadius, 
										iterations 	= 100);
	if(mlfit$error){ return(list(error=TRUE, errorMessage=mlfit$errorMessage)); }
	
	
	
	####################################
	# Calculate statistical significance 
	
	
	if(verbose) cat(sprintf("Calculating statistical significance under OUSS null hypothesis using Bernoulli trials..\n"));
	
	#Calculate statistical significance of absolute and relative power of detected peak with respect to fitted OUSS process
	significanceOUSS = significanceOfGlobalPeak(power_o		= mlfit$power_o, 
												lambda		= mlfit$lambda, 
												power_e		= mlfit$power_e, 
												time_step		= time_step, 
												frequencies	= frequencies[minPeakMode:length(frequencies)], 
												peakFreq	= peakFreq, 
												peakPower	= peakPower, 
												accuracy	= accuracy);
	
	#Correct estimated P-value
	significanceOUSS = interpolateWithinGrid(FAP_correction_grid, length(times), significanceOUSS, significanceOUSS);
	
	#Calculate significance of relative power of detected peak
	significanceOUSSlocal = significanceOfLocalPeak(power_o=mlfit$power_o, 
													lambda=mlfit$lambda, 
													power_e=mlfit$power_e, 
													time_step=time_step, 
													Nfreq=length(frequencies)-minPeakMode+1, 
													peakFreq=peakFreq, 
													peakPower=peakPower);
	
	#Return all results
	return(list(error			= FALSE,
				errorMessage	= "",
				frequencies		= frequencies, 
				periodogram 	= powers,
				fittedPS		= ps_ouss(frequencies, mlfit$power_o, mlfit$lambda, mlfit$power_e, time_step),
				peakMode		= peakMode,
				power_o			= mlfit$power_o, 
				lambda			= mlfit$lambda,
				power_e			= mlfit$power_e,
				time_step		= time_step,
				minPeakMode		= minPeakMode, 
				minFitMode		= minFitMode,
				MLL				= mlfit$MLL,
				P				= significanceOUSS,
				Plocal			= significanceOUSSlocal));
}
