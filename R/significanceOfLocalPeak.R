significanceOfLocalPeak <-
function(power_o, lambda, power_e, time_step, Nfreq, peakFreq, peakPower){
	Epower = ps_ouss(peakFreq, power_o, lambda, power_e, time_step);
	ratio = peakPower/Epower;
	if(is.nan(ratio)) return(NaN);
	return(1 - (1-exp(-ratio))^Nfreq);
}
