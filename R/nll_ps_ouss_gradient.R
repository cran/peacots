nll_ps_ouss_gradient <-
function(params, frequencies, powers, time_step){
	PSexpected 	= ps_ouss(frequencies, params[1]^2, abs(params[2]), params[3]^2, time_step);
	alpha 		= -1/PSexpected + powers/PSexpected^2;
	LL_grad		= c(0,0,0);
	LL_grad[1] 	= alpha %*% (2*params[1] 		* ps_ouss_gradient_wrt_power_o(frequencies, params[1]^2, abs(params[2]), params[3]^2, time_step));
	LL_grad[2]	= alpha %*% (sign(params[2]) 	* ps_ouss_gradient_wrt_lambda(frequencies, params[1]^2, abs(params[2]), params[3]^2, time_step));
	LL_grad[3]	= alpha %*% (2*params[3]		* ps_ouss_gradient_wrt_power_e(frequencies, params[1]^2, abs(params[2]), params[3]^2, time_step));
	return(-LL_grad);
}
