%% This function based on optimization of a
...model developed by Rakha et al 2012. The adj-Rsquared = 70%
...Rain intenisty must given in cm/h. However, this fucntion creates 
...a problem to the data as cm/h is a really high intensity of rain. Thus 
...we prefer to use a simple form of deceleration restricted to etab, 
...a parameter that tells us that the maximum deceleration 
...hampered by mechanical defeciencies of a vehicle per se.
function [bmx] = bmax(g, ri, mu, etab) 
    % bmx = g*(0.5088 - 0.03948 * RI); %(m/s)
    % bmx = g * etab * mu * (1.0 - 0.07759 * ri); %(m/s);
	 bmx = g * etab * mu;
end