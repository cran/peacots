nll_ps_ouss <-
function(params, frequencies, powers, time_step){
	PSexpected = ps_ouss(frequencies, params[1]^2, abs(params[2]), params[3]^2, time_step);
	LL = sum(- log(PSexpected) - powers/PSexpected);
	return(-LL);
}
