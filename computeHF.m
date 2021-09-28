function [veh] = computeHF(veh, h, xdist,ts)
% computeHF - Fund diagram task demand car following 
%
% USE 
% vehicle = computeHF(vehicle,h)
% in which vehicle is a struct with the following fields (default values
% given below)
%
% THIS function updates the parameters below:
% veh.TC      = 1;  % total task capacity
% veh.TD      = 0;  % total task demand
% veh.TS      = 0;  % task saturation
% veh.TDcf    = 0;  % car following task demand
% veh.TDdi    = 0;  % distraction task demand
% veh.SA      = 1;  % awareness/understanding level
% 
% % HF trajectory
% % [t,x,v,a] column vector 
% veh.trajHF = [veh.TDcf,veh.TDdi,veh.TS,veh.SA]; 
% 
% USING the FDTD functions specified below
%
% % default car following TDFD
% %         ^
% % tdmax |--
% %       |  :\
% %       |  : \
% %       |  :  \
% %   td0 |  :    -----
% %       |  :   :
% %        ------------> time headway
% %       0 hmin h0
% veh.FDTDcf.hmin  = 1.8; % s 
% veh.FDTDcf.tdmax = 1;   % -
% veh.FDTDcf.h0    = 2;   % s
% veh.FDTDcf.td0   = 0.4; % -
% 
% % 1st) Proposed Task Demand of Car Following
% %
% veh.FDTCcf.hmax  = 0.0; % [s] , maximum car following
% veh.FDTCcf.hcrit = 0.6; % [s] , critical car following 
% veh.FDTCcf.ha    = 2.0; % [s] , regular car following
% veh.FDTCcf.hb    = 4.0; % [s] , safe car follwoing
% veh.FDTDcf.tdmax = 1.0; % [-] , max task demand
% veh.FDTDcf.tdcrit= 0.8; % [-] , critical task demand
% veh.FDTDcf.tda   = 0.6; % [-] , safe task demand
% veh.FDTDcf.tdb   = 0.4; % [-] , free flow task demand
%
% % 2nd) Proposed Task Demand of Car Following
% veh.FDTCcf.hexp  = 10.00; % [s] , maximum car following
% 
% % default distraction TDFD
% %       ^
% %       |
% % tdmax |..........___    
% %       |         /: :\
% %       |       /  : : \
% %       |     /    : :  \
% %       |   /      : :   \
% %        ---------------------> space
% %         x0       0 x1  x2
% veh.FDTDdi.tdmax  = 0.8; % s
% veh.FDTDdi.x0     = -400; % m distance to distraction from which driver starts being distracted 
% veh.FDTDdi.x1     = 200;  % m distance in full view of distraction
% veh.FDTDdi.x2     = 400;  % m distance after distraction in which driver fully recovers
% 
% % Awareness relations
% % 
% %         ^
% %   samax |------
% %         |     : \
% %         |     :  \
% %         |     :   \
% %   samin |.......... ----
% %         |     :    :
% %          -------------> task saturation
% %         0  tscrit tsmax
% 
% % function that governs perception errors
% veh.SAFn.samax = 1;
% veh.SAFn.samin   = 0.5;
% veh.SAFn.tscrit  = 0.95;
% veh.SAFn.tsmax   = 2;
% % reaction time effect?
% veh.SAFn.taumax  = 0.4; % seconds
%
%Outputs
%   vehicle    updated vehicle structure

% if no distraction task
if nargin<3
    xdist = inf;
end;

% if ts is given as input, we do not relate SA to TD, but simply compute
% the result. In most cases, the three inputs are now vectors so that the
% result can be used to make a graph
demoMode = nargin == 4;

%% First driving task
if ~demoMode && isinf(h)
    h = 60; %s
end;

% % we determine the hmin parameter in the CF FDTD, which we increase in case
% % of large accelerations:
%  if veh.a<-veh.par.bcomf % acc < comfortable deceleration?
%      % we linearly increase it
%      hminfac = (veh.a+veh.par.bcomf)./(veh.par.bmax-veh.par.bcomf);
%  else
%      hminfac = 0;
%  end
 hminfac = 0;
 hmin  = (1+hminfac).*veh.FDTDcf.hmin;   % s 

% we scale the FD according task capacity
 tdmax = veh.FDTDcf.tdmax  ./ veh.TC;   % -
 h0    = veh.FDTDcf.h0;                 % s
 td0   = veh.FDTDcf.td0  ./ veh.TC;     % -

% we also scale the FDTD according to additional task difficulty due to Rain Intensity
 tdmax = tdmax .* veh.FDTDcf.D;   % -
 td0   = td0   .* veh.FDTDcf.D;   % -
 
% 1st proposal
% hmax  = veh.FDTCcf.hmax;  % [s] , critical car following
% hcrit = veh.FDTCcf.hcrit;  % [s] , critical car following 
% ha    = veh.FDTCcf.ha;     % [s] , regular car following
% hb    = veh.FDTCcf.hb;     % [s] , safe car follwoing
% tdmax  = veh.FDTDcf.tdmax ./ veh.TC;  % [-] , max task demand
% tdcrit = veh.FDTDcf.tdcrit ./ veh.TC; % [-] , critical task demand
% tda    = veh.FDTDcf.tda ./ veh.TC;   % [-] , safe task demand
% tdb    = veh.FDTDcf.tdb ./ veh.TC;   % [-] , free flow task demand

% 2nd proposal
% hexp  = veh.FDTCcf.hexp;  % [s]

% Task demand equals
 veh.TDcf = (h<=hmin) .* tdmax + ...
     (h>hmin & h<=h0) .* (tdmax - (h-hmin)*(tdmax-td0)./(h0-hmin)) + ...
               (h>h0) .* td0 ;

% Task Demand Equals
% veh.TDcf = (h<=hcrit) .* (tdmax - ((h-hcrit)/(hmax - hcrit)) * (tdmax - tdcrit)) * veh.FDTDcf.D + ...
%     (h>hcrit & h<=ha) .* (tdcrit - ((h-ha)/(hcrit - ha)) * (tdcrit - tda)) * veh.FDTDcf.D + ...
%         (h>ha & h<hb) .* (tda - ((h-hb)/(ha - hb))*(tda - tdb)) * veh.FDTDcf.D +...
%                (h>hb) .* tdb * veh.FDTDcf.D;

% Task Demand Equals
% veh.TDcf = exp(-h/hexp) * veh.FDTDcf.D;

%% distraction task
if ~demoMode && isinf(xdist)
    xdist = 10e3;
end

% shortcuts
% we scale the FD according task capacity
tdmax = veh.FDTDdi.tdmax / veh.TC; 
x0 = veh.FDTDdi.x0; 
x1 = veh.FDTDdi.x1; 
x2 = veh.FDTDdi.x2; 
x  = xdist;

% negative distance means we're are closing in, positive distance means
% we're moving away from the distraction
veh.TDdi = 0 + (x > x0 & x<= 0) .* (1-(x./x0)) .* tdmax + ... % both x and x0 are negative
               (x >  0 & x<=x1) .* tdmax + ... % maximum distraction        
               (x > x1 & x<=x2) .* (1-(x-x1)./(x2-x1)) .* tdmax;

%% total task demand & task saturation
if ~demoMode
    veh.TD = veh.TDcf + veh.TDdi;
    veh.TS = veh.TD ./ veh.TC;
else
    veh.TD = ts;
    veh.TS = ts;
end;

%% SA / understanding & reaction time
% shortcuts
samax   = veh.SAFn.samax;
samin   = veh.SAFn.samin;
tscrit  = veh.SAFn.tscrit;
tsmax   = veh.SAFn.tsmax;
ts      = veh.TS;

% resulting sa level
veh.SA =     (ts<=tscrit) + ...
  (ts>tscrit & ts<=tsmax) .* (samax - (ts-tscrit)./(tsmax-tscrit).*(samax-samin)) + ...
               (ts>tsmax) .* samin;    
