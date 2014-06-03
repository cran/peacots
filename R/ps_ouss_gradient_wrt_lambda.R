ps_ouss_gradient_wrt_lambda <-
function(freq, power_o, lambda, power_e, time_step){
	rho = exp(-lambda*time_step);
	return(2*power_o*rho*time_step*(1-rho)*(rho^2+rho+cos(2*pi*freq*time_step)*(1-3*rho))/(1+rho^2-2*rho*cos(2*pi*freq*time_step))^2);
}
