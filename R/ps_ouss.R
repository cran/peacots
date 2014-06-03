ps_ouss <-
function(freq, power_o, lambda, power_e, time_step){
	rho = exp(-lambda*time_step);
	return(power_o*(1-rho)^2/(1+rho^2-2*rho*cos(2*pi*freq*time_step)) + power_e);
}
