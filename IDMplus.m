function [acc, acc_nobound] = IDMplus(vehicle, ds, dv, amax, v0, s0, T, bcomf, alfa, bmax)
% IDMPLUS - car following model according to Schakel & van Arem (2010)
%
% USE It can be called in two ways
% (1) acc = IDMplus(v, ds, dv, amax, v0, s0, T, bcomf, alfa, bmax)
%     in which all inputs are 1xN sized vectors (depicting N vehicles)
%     v: vehicle speedsw, ds: distance gaps, dv: speed differences
%
% (2) acc = IDMplus(vehicle, ds, dv)
%     in which vehicle is a struct with the following fields
%     vehicle.v                  % prevailing speed
%     vehicle.par.amax = 3;      % m/s2 max acc
%     vehicle.par.v0   = 35;     % m/s  des speed
%     vehicle.par.s0   = 5;      % m    stopping dist
%     vehicle.par.T    = 1.2;    % s    min time headway
%     vehicle.par.bcomf= 3;      % m/s2 comfortable acceleration
%     vehicle.par.alfa = 4;      % -    scale parameter (power 1st term)
%     vehicle.par.bmax = -8;     % m/s2 max deceleration
%
%     NB: v0 and T are multiplied by
%     vehicle.ctl.v0fac
%     vehicle.ctl.Tfac
% 
% (3) In order to fully comply with the rationale of IDM,IDM+
%     we introduce the following equation in which; The acceleration takes
%     into account the mechanical capabilities of a vehicle and the
%     preferences of a driver for acceleration. The minimization term combines both preferences and mechanecal capabilities. 
% 
%    amax = min(amax_mech, amax_beh)
% 
% Outputs
%     acc ; vehicle acceleration according to IDM+

%%

% parameters per vehicle
if isstruct(vehicle)
    v    = vehicle.v;
    amax_beh = vehicle.par.amax_beh;                 % m/s2 max acc
    amax_mech = vehicle.par.amax_mech;               % m/s2 max acc
    v0   = vehicle.ctl.v0fac * vehicle.par.v0;       % m/s  des speed
    s0   = vehicle.par.s0;                           % m    stopping dist
    T    = vehicle.ctl.Tfac * vehicle.par.T;         % s    min time headway
    bcomf= vehicle.par.bcomf;                        % m/s2 comfortable acceleration
    alfa = vehicle.par.alfa;                         % -    scale parameter (power 1st term)
    bmax = vehicle.par.bmax;                         % m/s2 max deceleration, 8 m/s2 it the value here
%    cr   = vehicle.par.Cr;
%    c2   = vehicle.par.c2;
%    c3   = vehicle.par.c3;
%    g    = vehicle.par.g;
%    mu   = vehicle.par.mu;
%    percm= vehicle.par.percm;
elseif nargin == 9
    v    = vehicle;
else
    error('IDMplus: incorrect nr of inputs');
end

amax = min(amax_mech, amax_beh);

% IDM+
ds_star= s0 + max(0, T.*v + v.*dv ./ (2.*sqrt(amax.*bcomf)));
acc_nobound = amax * ( min( 1-(v./v0).^alfa, 1-(ds_star./ds).^2) );
acc = max( acc_nobound, -bmax );


