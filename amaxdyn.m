%% By this function a dynamic calculation of the deceleration of a vehicle 
...is derived. The parameters and constants of the functions are the followings
...
function [amxd] = amaxdyn(Cr, c3, g, Fmax, m)
   amxd = (Fmax - g * m * Cr * c3 *10^(-3)) / m;
end
