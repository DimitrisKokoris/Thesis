function [ veh ] = createDefVehicle(Sim)
%createDefVehicle: create default vehicle with default driving and HF
%                  params
%
%   Detailed explanation goes here

% fields / flags for simulation
veh.ID        = 1;      % ID so the vehicle can be tagged
veh.tgendes   = NaN;    % time instant vehicle should be generated   
veh.kgendes   = NaN;    % index time instant vehicle should be generated   
veh.tgen      = NaN;    % time instant vehicle was actually generated
veh.kgen      = NaN;    % index time instant vehicle was actually generated
veh.kfin      = NaN;    % index time instant vehicle finished his course
veh.tqueue    = 0;      % amount of time vehicle spent while waiting to get on the track 
veh.active    = false;  % flag to mark vehicle is/was active (i.e. driving)
veh.finished  = false;  % flag to mark vehicle has left the network

% graphics
veh.hobj      = -1;     % handle to graphics object(s)
veh.color     = 'k';    % color trajectory;
veh.width     = 0.1;    % thickness trajectory;

% incidents
veh.incident  = false;  % flags if the vehicle was part of an incident
veh.inc.hobj  = -1;     % handle to graphics object
veh.inc.x     = 0;      % location
veh.inc.t     = 0;      % time

% default vehicle params IDM(+)
veh.par.tau  = 0;      % s    reaction time
veh.par.amax_beh = 3;      % m/s2 max acc
veh.par.v0   = 33.33;  % m/s  des speed
veh.par.s0   = 8;      % m    stopping dist
veh.par.T    = 1.2;    % s    min time headway
veh.par.bcomf= 2;      % m/s2 comfortable acceleration
veh.par.alfa = 4;      % -    scale parameter (power 1st term)
%veh.par.bmax = 8;           % m/s2 max deceleration
veh.par.ant  = 1;      % nr of vehicles used for anticipation NOT IMPLEMENTED

% proposed method to account for environmental conditions(rain)
veh.par.mu    = Sim.mu;  % - coefficient of friction
veh.par.ri    = Sim.ri;  % mm/h rain intensity per hour
veh.par.g     = 9.80665; % m/s^2 gravitational acceleration
veh.par.etab  = 1.0;     % - braking efficiency
veh.par.percm = 0.55;     % kg mass of a vehicle on the tracive axle
veh.par.Cr    = 1.25;    % - rolling coefficient of good asphalt condition
	% Rolling resistant coefficient
veh.par.c2    = 0.0328;  % - Tire type
veh.par.c3    = 4.575;   % - Tire type
	% Engine characteristics and vehicle mass
veh.par.P     = 90;      % kW power of a vehicle in
veh.par.eta   = 0.94;    % - engine efficiency, typical values are between 0.89 and 0.94
veh.par.m     = 1400;    % kg mass of a vehicle expressed
	
% vehicle response control factors
veh.ctl.v0fac    = 1; % multiplicator for free speed
veh.ctl.Tfac     = 1; % multiplicator for desired headway
veh.ctl.tauextra = 0; % additional reaction time
veh.ctl.antfac   = 1; % addition OR reduction in anticipation (minimum = 1) NOT IMPLEMENTED

% vehicle kinematics
veh.t = 0;
veh.x = 0;
veh.v = veh.par.v0;
veh.a = 0;

% GROSS gap and headway
veh.s = 0;
veh.h = 0;

% Compute forces
% Ft   = Ftract(veh.par.eta, veh.par.P, veh.v);
Fmx = Fmax(veh.par.g, veh.par.m, Sim.mu, veh.par.percm);

% Here call the function of bmax, since the bmax is not dependent on the speed, only on RI though, we call here the bmax function
veh.par.bmax= bmax(veh.par.g, veh.par.ri,veh.par.mu, veh.par.etab);       % Calculate only once.
veh.par.amax_mech= amaxdyn(veh.par.Cr, veh.par.c3, veh.par.g, Fmx, veh.par.m); % Calculate only once.

% vehicle trajectory
% [t,x,v,a,s,h] column vector 
veh.traj = [veh.t,veh.x,veh.v,veh.a,0,0]; 

% Human factors
veh.TC      = 1.0;  % total task capacity
veh.TD      = 0;    % total task demand
veh.TS      = 0;    % task saturation
veh.TDcf    = 0;    % car following task demand
veh.TDdi    = 0;    % distraction task demand
veh.SA      = 1;    % awareness/understanding level
% vehicle perception biases
veh.percbias = -1; % [-1, 1]: [underestimation, overestimation]

% HF trajectory
% [t,x,v,a] column vector 
veh.trajHF = [veh.t,veh.x,veh.TDcf,veh.TDdi,veh.SA,veh.par.tau,NaN]; 

% Default car following TDFD
%         ^
% tdmax |--
%       |  :\
%       |  : \
%       |  :  \
%   td0 |  :    -----
%       |  :   :
%        ------------> time headway
%       0 hmin h0

if Sim.doHF_TD == true
    veh.FDTDcf.f     = 2.0/1e2;                                 % Capacity reduction factor as Calvert et al 2013 proposed DOI: 10.1109/ITSC.2013.6728439
    veh.FDTDcf.D     = round(1/(1 - veh.FDTDcf.f * Sim.ri ),2); % Task Difficulty due to rain! D = 1/h_N / 1/h_R = h_R / h_N = 1/(1 - veh.FDTDcf.f * Sim.ri )
else
    veh.FDTDcf.f     = 0;                        
    veh.FDTDcf.D     = round(1/(1 - veh.FDTDcf.f * Sim.ri ),2); 
end
veh.FDTDcf.hmin  = 0.8;                            % s 
veh.FDTDcf.tdmax = 1;                              % -
veh.FDTDcf.h0    = 3;                              % s
veh.FDTDcf.td0   = 0.5;                            % -

% veh.FDTCcf.hmax  = 0.0; % [s] , maximum car following
% veh.FDTCcf.hcrit = 0.6; % [s] , critical car following 
% veh.FDTCcf.ha    = 2.0; % [s] , regular car following
% veh.FDTCcf.hb    = 4.0; % [s] , safe car follwoing
% veh.FDTDcf.tdmax = 1.0; % [-] , max task demand
% veh.FDTDcf.tdcrit= 0.8; % [-] , critical task demand
% veh.FDTDcf.tda   = 0.6; % [-] , safe task demand
% veh.FDTDcf.tdb   = 0.4; % [-] , free flow task demand

% veh.FDTCcf.hexp    = 8;   

% default distraction TDFD
%       ^
%       |
% tdmax |..........___    
%       |         /: :\
%       |       /  : : \
%       |     /    : :  \
%       |   /      : :   \
%        ---------------------> space
%         x0       0 x1  x2
veh.FDTDdi.tdmax  = 0.8; % s
veh.FDTDdi.x0     = -400; % m distance to distraction from which driver starts being distracted 
veh.FDTDdi.x1     = 200;  % m distance in full view of distraction
veh.FDTDdi.x2     = 400;  % m distance after distraction in which driver fully recovers

% Awareness relations
% 
%         ^
%   samax |------
%         |     : \
%         |     :  \
%         |     :   \
%   samin |.......... ----
%         |     :    :
%          -------------> task saturation
%         0  tscrit tsmax

% function that governs perception errors
veh.SAFn.samax   = 1;
veh.SAFn.samin   = 0.5;
veh.SAFn.tscrit  = 0.8;
veh.SAFn.tsmax   = 2;
% reaction time effect? This is a composition of sensation time and
% perception time tau = tau_{senstation} + tau_{perception}
veh.SAFn.taumax  = 0.5; % seconds. 

end

