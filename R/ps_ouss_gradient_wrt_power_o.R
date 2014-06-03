ps_ouss_gradient_wrt_power_o <-
function(freq, power_o, lambda, power_e, time_step){
	rho = exp(-lambda*time_step);
	return((1-rho)^2/(1+rho^2-2*rho*cos(2*pi*freq*time_step)));
}
