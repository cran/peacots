LombScarglePeriodogram <-
function(times, signal){
	N = length(times);
	if(N<2) return(list(freq=c(), spec=c()));
	T = times[N] - times[1];
	
	#Get average, minimum and maximum time steps and see if time increments are uniform
	meandt=T/(N-1);
	mindt=maxdt=times[2]-times[1];
	for(n in 1:(N-1)){
		mindt = min(mindt, times[n+1]-times[n]);
		maxdt = max(maxdt, times[n+1]-times[n]);
	}
	uniform = (maxdt-mindt < 0.00001*meandt);
	
	f0 = 1.0/(N*meandt);
	maxMode = floor(N/2); #upper limit given by Nyquist frequency
	spec = freq = (0:maxMode)*f0;
	
	if(uniform){
		#time steps are uniform, so use classical DFT formula for optimization reasons
		time_indices = 0:(N-1);
		for(m in 0:maxMode){
			FRe = mean(signal* cos(-(2.0*pi*m*time_indices)/N));
			FIm = mean(signal* sin(-(2.0*pi*m*time_indices)/N));
			spec[m+1] = T*(FRe*FRe + FIm*FIm);
		}
		
	}else{
		#time steps not uniform, so use more complicated Lomb-Scargle periodogram formula
		Xmean = mean(signal);
		spec[1] = Xmean*Xmean*T;
		for(m in 1:maxMode){
			f = freq[m+1];
			
			#calculate time shift
			C = sum(cos(4*pi*f*times));
			S = sum(sin(4*pi*f*times));
			tau = atan(S/C)/(4*pi*f);
			
			#calculate power
			C = sum((cos(2*pi*f*(times-tau)))^2);
			S = sum((sin(2*pi*f*(times-tau)))^2);
			FC = sum((signal-Xmean)*cos(-2*pi*f*(times-tau)))/sqrt(2*C*N);
			FS = sum((signal-Xmean)*sin(-2*pi*f*(times-tau)))/sqrt(2*S*N);
			spec[m+1] = T*(FC*FC + FS*FS);
		}
	}
	
	return(list(freq=freq, spec=spec));
}
