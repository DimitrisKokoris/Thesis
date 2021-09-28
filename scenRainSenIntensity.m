%% General sim preparations for Rain Conditions
% clear up
clear all variables; clc;

dofig  = false;
dosavedata = true;
dosavefig  = true;

% timeline
Sim.dt       = 0.1;             % sec
Sim.Tsim     = 15*60;           % sec
Sim.timeline = 0:Sim.dt:Sim.Tsim-Sim.dt;    % sec
Sim.numsteps = numel(Sim.timeline);
% update progress every ... steps
Sim.dtprogress = 450;
% road stretch
Sim.xstart =  0;
Sim.xmax   = 6e3;  % exit at 3 km
% max nr of vehicles generated:
Sim.numvehicles = 350;
% Safe zone: we cannot apply reaction time logic close to origin.
Sim.safezone = 5*40; % 5 (sec) * 40 (m/s) (ultra long reaction time * very high speed)

%% Scenario preparations

% base OD profile (sine curve)
q = zeros(size(Sim.timeline));
sinebase  = 1+sin(1.5*pi/Sim.Tsim * Sim.timeline);
% qbase = 900 * sinebase;
qbase = [900  * ones(1,0.1*Sim.numsteps), ...
    2200 * ones(1,0.3*Sim.numsteps), ...
    900  * ones(1,0.6*Sim.numsteps)];

%qbase = 1400  * ones(1,Sim.numsteps);
%q = poissrnd(qbase);
Sim.q = qbase;
%% here's the list of scenario's
scenarios = { ...
% BASE CASE
% 	     'scenBase_I=0.0',...
% BASE CASE AND MVD
%   	   'scenBaseR_MVD_I=0.0',...
%        'scenBaseR_MVD_I=0.5',...
%        'scenBaseR_MVD_I=1.0',...
%        'scenBaseR_MVD_I=2.0',...
%        'scenBaseR_MVD_I=3.0',...
%        'scenBaseR_MVD_I=4.0',...
%        'scenBaseR_MVD_I=5.0',...
%        'scenBaseR_MVD_I=6.0',...
%        'scenBaseR_MVD_I=7.0',...
%        'scenBaseR_MVD_I=8.0',...
%        'scenBaseR_MVD_I=9.0',...
% 	   'scenBaseR_MVD_I=10.0',...
% BASE CASE + MVD + INDUCE TASK DIFFICULTY + REACTION TIME
%        'scenBaseR_MVD_TD_RT_I=0.0',...',...
%        'scenBaseR_MVD_TD_RT_I=0.5',...
%        'scenBaseR_MVD_TD_RT_I=1.0',...
%        'scenBaseR_MVD_TD_RT_I=2.0',...
%        'scenBaseR_MVD_TD_RT_I=3.0',...
%        'scenBaseR_MVD_TD_RT_I=4.0',...
%        'scenBaseR_MVD_TD_RT_I=5.0',...
%        'scenBaseR_MVD_TD_RT_I=6.0',...
%        'scenBaseR_MVD_TD_RT_I=7.0',...
%        'scenBaseR_MVD_TD_RT_I=8.0',...
%        'scenBaseR_MVD_TD_RT_I=9.0',...
% 	   'scenBaseR_MVD_TD_RT_I=10.0',...  
% Scenarios with 75% underestimation and 25% overestimation
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s) 
%        'scenBaseR_MVD_TD_PES_I=0.0',...
%        'scenBaseR_MVD_TD_PES_I=0.5',...
%        'scenBaseR_MVD_TD_PES_I=1.0',...
%        'scenBaseR_MVD_TD_PES_I=2.0',...
%        'scenBaseR_MVD_TD_PES_I=3.0',...
%        'scenBaseR_MVD_TD_PES_I=4.0',...
%        'scenBaseR_MVD_TD_PES_I=5.0',...
%        'scenBaseR_MVD_TD_PES_I=6.0',...
%        'scenBaseR_MVD_TD_PES_I=7.0',...
%        'scenBaseR_MVD_TD_PES_I=8.0',...
%        'scenBaseR_MVD_TD_PES_I=9.0',...
% 	   'scenBaseR_MVD_TD_PES_I=10.0',...
% Scenarios with 75% underestimation and 25% overestimation
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(Dv) 
%        'scenBaseR_MVD_TD_PEV_I=0.0',...
%        'scenBaseR_MVD_TD_PEV_I=0.5',...
%        'scenBaseR_MVD_TD_PEV_I=1.0',...
%        'scenBaseR_MVD_TD_PEV_I=2.0',...
%        'scenBaseR_MVD_TD_PEV_I=3.0',...
%        'scenBaseR_MVD_TD_PEV_I=4.0',...
%        'scenBaseR_MVD_TD_PEV_I=5.0',...
%        'scenBaseR_MVD_TD_PEV_I=6.0',...
%        'scenBaseR_MVD_TD_PEV_I=7.0',...
%        'scenBaseR_MVD_TD_PEV_I=8.0',...
%        'scenBaseR_MVD_TD_PEV_I=9.0',...
% 	     'scenBaseR_MVD_TD_PEV_I=10.0',...
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s) 
%        'scenBaseR_MVD_TD_PESO_I=0.0',...
%        'scenBaseR_MVD_TD_PESO_I=0.5',...
%        'scenBaseR_MVD_TD_PESO_I=1.0',...
%        'scenBaseR_MVD_TD_PESO_I=2.0',...
%        'scenBaseR_MVD_TD_PESO_I=3.0',...
%        'scenBaseR_MVD_TD_PESO_I=4.0',...
%        'scenBaseR_MVD_TD_PESO_I=5.0',...
%        'scenBaseR_MVD_TD_PESO_I=6.0',...
%        'scenBaseR_MVD_TD_PESO_I=7.0',...
%        'scenBaseR_MVD_TD_PESO_I=8.0',...
%        'scenBaseR_MVD_TD_PESO_I=9.0',...
% 	   'scenBaseR_MVD_TD_PESO_I=10.0',...
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(Dv) 
%        'scenBaseR_MVD_TD_PEVO_I=0.0',...
%        'scenBaseR_MVD_TD_PEVO_I=0.5',...
%        'scenBaseR_MVD_TD_PEVO_I=1.0',...
%        'scenBaseR_MVD_TD_PEVO_I=2.0',...
%        'scenBaseR_MVD_TD_PEVO_I=3.0',...
%        'scenBaseR_MVD_TD_PEVO_I=4.0',...
%        'scenBaseR_MVD_TD_PEVO_I=5.0',...
%        'scenBaseR_MVD_TD_PEVO_I=6.0',...
%        'scenBaseR_MVD_TD_PEVO_I=7.0',...
%        'scenBaseR_MVD_TD_PEVO_I=8.0',...
%        'scenBaseR_MVD_TD_PEVO_I=9.0',...
% 	   'scenBaseR_MVD_TD_PEVO_I=10.0',...
% Scenarios with 75% underestimation and 25% overestimation
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) 
%        'scenBaseR_MVD_TD_PESV_I=0.0',...
%        'scenBaseR_MVD_TD_PESV_I=0.5',...
%        'scenBaseR_MVD_TD_PESV_I=1.0',...
%        'scenBaseR_MVD_TD_PESV_I=2.0',...
%        'scenBaseR_MVD_TD_PESV_I=3.0',...
%        'scenBaseR_MVD_TD_PESV_I=4.0',...
%        'scenBaseR_MVD_TD_PESV_I=5.0',...
%        'scenBaseR_MVD_TD_PESV_I=6.0',...
%        'scenBaseR_MVD_TD_PESV_I=7.0',...
%        'scenBaseR_MVD_TD_PESV_I=8.0',...
%        'scenBaseR_MVD_TD_PESV_I=9.0',...
% 	     'scenBaseR_MVD_TD_PESV_I=10.0',...
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% % BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) 
%        'scenBaseR_MVD_TD_PESVO_I=0.0',...
%        'scenBaseR_MVD_TD_PESVO_I=0.5',...
%        'scenBaseR_MVD_TD_PESVO_I=1.0',...
%        'scenBaseR_MVD_TD_PESVO_I=2.0',...
%        'scenBaseR_MVD_TD_PESVO_I=3.0',...
%        'scenBaseR_MVD_TD_PESVO_I=4.0',...
%        'scenBaseR_MVD_TD_PESVO_I=5.0',...
%        'scenBaseR_MVD_TD_PESVO_I=6.0',...
%        'scenBaseR_MVD_TD_PESVO_I=7.0',...
%        'scenBaseR_MVD_TD_PESVO_I=8.0',...
%        'scenBaseR_MVD_TD_PESVO_I=9.0',...
	     'scenBaseR_MVD_TD_PESVO_I=10.0',...
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(v)
%        'scenBaseR_MVD_TD_PESV_BAS_I=0.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=0.5',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=1.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=2.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=3.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=4.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=5.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=6.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=7.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=8.0',...
%        'scenBaseR_MVD_TD_PESV_BAS_I=9.0',...
% 	    'scenBaseR_MVD_TD_PESV_BAS_I=10.0',... 
 % BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(T)
%        'scenBaseR_MVD_TD_PESV_BAT_I=0.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=0.5',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=1.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=2.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=3.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=4.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=5.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=6.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=7.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=8.0',...
%        'scenBaseR_MVD_TD_PESV_BAT_I=9.0',...
% 	     'scenBaseR_MVD_TD_PESV_BAT_I=10.0'...  
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(v,T)
%         'scenBaseR_MVD_TD_PESV_BAST_I=0.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=0.5',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=1.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=2.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=3.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=4.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=5.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=6.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=7.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=8.0',...
%         'scenBaseR_MVD_TD_PESV_BAST_I=9.0',...
% 	    'scenBaseR_MVD_TD_PESV_BAST_I=10.0'...  
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +
% BEHAVIORAL ADAPTATION(v,T) + HETEROGENEITY
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=0.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=0.5',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=1.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=2.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=3.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=4.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=5.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=6.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=7.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=8.0',...
%        'scenBaseR_MVD_TD_PESV_BAST_H_I=9.0',...
%  	   'scenBaseR_MVD_TD_PESV_BAST_H_I=10.0'... 
};

% run over all of them
numscen = numel(scenarios);

%  I must fix this!!!(rection time)
for itau = 0
    Sim.base_tau    = itau/10;
    tau_str = sprintf('tau=%.2d',itau);
    
    fprintf('\n-------------------------\n');
    fprintf('BASE REACTION TIME = %.1f  \n',Sim.base_tau);
    fprintf('-------------------------\n');
    
    % over all the scenarios
    for iscenario = 1 : numscen
        
        % BASE Settings for the scenarios
        
       % A bottleneck at 4 and 4.5
        Sim.xDisturbance  = [4e3,4.5e3];    % disturbance location
        Sim.tDisturbance  = [1, 15] * 60;   % disturbance
        Sim.disturbanceStrength = 0.36;      % Severity (qcap = 2100 veh/h, u = 33.33 m/s)
    
         % Base scenario:
        Sim.doHF        = true;      % do we compute HF?
        Sim.doDistraction = false;   % do we simulate a distraction?
        Sim.doDisturbance = true;    % do we force an additional disturbances, lets assume a flow conserving bottleneck?
        Sim.doHeterogeneity = false; % do we vary with critical task saturation?
        %Sim.base_tau    = itau;     % do we give drivers a base reaction time?
        Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
        Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
        Sim.doHF_tau    = false;     % do we have tau = tau_{senstation} + tau_{perception}
        Sim.doHF_v0     = false;     % do we have HF affect v0
        Sim.doHF_T      = false;     % do we have HF affect T
        Sim.doHF_TD     = false;     % do we have HF affect Task Difficulty
        
        % What is the perception biase of the driving population?
        Sim.doHF_prob_overestimation = 0.25;
        % selected vehicles for HF plots
        Sim.iHF = [191,192];
        
        % select scenario
        Sim.scenario = scenarios{iscenario};
        switch Sim.scenario
		%% BASE CASE SCENARIO
		   case 'scenBase_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
        %% MVD SET OF SCENARIOS
		   case 'scenBaseR_MVD_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_tau    = true; % This regrads the base reaction time of the drivers
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
               
           case 'scenBaseR_MVD_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
                
           case 'scenBaseR_MVD_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);

    %% MECHANICAL VEHICLE DYNAMICS + INDUCED TASK DIFFICULTY + REACTION TIME
		   case 'scenBaseR_MVD_TD_RT_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_RT_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_RT_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_RT_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_RT_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_RT_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_RT_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_RT_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_RT_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_RT_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_RT_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_RT_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_tau    = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);

    %% MECHANICAL VEHICLE DYNAMICS + INDUCED TASK DIFFICULTY + PERCEPTION ERRORS DISTANCE
		   case 'scenBaseR_MVD_TD_PES_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PES_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PES_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PES_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PES_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PES_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PES_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_PES_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_PES_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_PES_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PES_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_PES_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);

%% MECHANICAL VEHICLE DYNAMICS + INDUCED TASK DIFFICULTY + PERCEPTION ERRORS SPEED
		   case 'scenBaseR_MVD_TD_PEV_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEV_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                SSim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEV_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEV_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                    Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                    Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEV_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEV_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEV_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_PEV_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_PEV_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_PEV_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEV_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_PEV_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s) 
		   case 'scenBaseR_MVD_TD_PESO_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESO_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESO_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESO_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESO_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESO_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESO_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_PESO_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_PESO_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_PESO_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESO_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_PESO_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = false;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(Dv) 

		   case 'scenBaseR_MVD_TD_PEVO_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEVO_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                SSim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEVO_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEVO_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                    Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                    Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEVO_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEVO_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PEVO_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_PEVO_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_PEVO_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_PEVO_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PEVO_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_PEVO_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = false;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenarios with 75% underestimation and 25% overestimation
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) 
            case 'scenBaseR_MVD_TD_PESV_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESV_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenarios with 75% OVERESTIMATION and 25% UNDERESTIMATION
% BASE CASE + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) 
            case 'scenBaseR_MVD_TD_PESVO_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
            case 'scenBaseR_MVD_TD_PESVO_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_prob_overestimation = 0.75;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(v)
          case 'scenBaseR_MVD_TD_PESV_BAS_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAS_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAS_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAS_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);    
           case 'scenBaseR_MVD_TD_PESV_BAS_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAS_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAS_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAS_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAS_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);   
           case 'scenBaseR_MVD_TD_PESV_BAS_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAS_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAS_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(T)
		            case 'scenBaseR_MVD_TD_PESV_BAT_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAT_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAT_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAT_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);    
           case 'scenBaseR_MVD_TD_PESV_BAT_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T    = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAT_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAT_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAT_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0  
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAT_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);   
           case 'scenBaseR_MVD_TD_PESV_BAT_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAT_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAT_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_T     = true;    % do we have HF affect v0 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +  BEHAVIORAL ADAPTATION(v,T)
          case 'scenBaseR_MVD_TD_PESV_BAST_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAST_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAST_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);
           case 'scenBaseR_MVD_TD_PESV_BAST_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true; 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);    
           case 'scenBaseR_MVD_TD_PESV_BAST_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true; 
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAST_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAST_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAST_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0  
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
           case 'scenBaseR_MVD_TD_PESV_BAST_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);   
           case 'scenBaseR_MVD_TD_PESV_BAST_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAST_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD);  
           case 'scenBaseR_MVD_TD_PESV_BAST_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;    % do we have HF affect reaction time (attention lag)
                Sim.doHF_v0     = true;    % do we have HF affect v0 
                Sim.doHF_T     = true;
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity,...
                Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0,...
                Sim.doHF_T, Sim.doHF_TD); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASE AND BASE REACTION TIME, TASK DIFFICULTY, PERCEPTION ERRORS, HETEROGENEITY, BEHAVIORAL ADAPTATION SPEED AND HEADWAY!
		   case 'scenBaseR_TD_PET_H_BAST_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_TD_PET_H_BAST_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_TD_PET_H_BAST_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_TD_PET_H_BAST_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_TD_PET_H_BAST_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;               
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_TD_PET_H_BAST_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_TD_PET_H_BAST_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;               % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;               % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;               % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;           % do we vary with critical task saturation?
                Sim.doHF_v0     = true;               % do we have HF affect v0
                Sim.doHF_T      = true;               % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_TD_PET_H_BAST_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_TD_PET_H_BAST_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_TD_PET_H_BAST_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_TD_PET_H_BAST_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;               % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;               % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;               % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;           % do we vary with critical task saturation?
                Sim.doHF_v0     = true;               % do we have HF affect v0
                Sim.doHF_T      = true;               % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_TD_PET_H_BAST_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;      % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;      % do we have HF affect perception of speed differences
                Sim.doHF_tau    = false;      % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASE CASE  + MVD + INDUCE TASK DIFFICULTY +  PERCEPTION ERRORS(s,Dv) +
% BEHAVIORAL ADAPTATION(v,T) + HETEROGENEITY
		   case 'scenBaseR_MVD_TD_PESV_BAST_H_I=0.0'
                Sim.mu = 0.60;
                Sim.ri = 0.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESV_BAST_H_I=0.5'
                Sim.mu = 0.40;
                Sim.ri = 0.5;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=1.0'
                Sim.mu = 0.40;
                Sim.ri = 1.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=2.0'
                Sim.mu = 0.40;
                Sim.ri = 2.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=3.0'
                Sim.mu = 0.40;
                Sim.ri = 3.0;               
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESV_BAST_H_I=4.0'
                Sim.mu = 0.40;
                Sim.ri = 4.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=5.0'
                Sim.mu = 0.40;
                Sim.ri = 5.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;               % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;               % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;               % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;           % do we vary with critical task saturation?
                Sim.doHF_v0     = true;               % do we have HF affect v0
                Sim.doHF_T      = true;               % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
               
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=6.0'
                Sim.mu = 0.40;
                Sim.ri = 6.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=7.0'
                Sim.mu = 0.40;
                Sim.ri = 7.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
                
           case 'scenBaseR_MVD_TD_PESV_BAST_H_I=8.0'
                Sim.mu = 0.40;
                Sim.ri = 8.0;
                Sim.doHF_TD     = true;
                Sim.doHF_perc_s = true;     % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;     % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;     % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
           
		   case 'scenBaseR_MVD_TD_PESV_BAST_H_I=9.0'
                Sim.mu = 0.40;
                Sim.ri = 9.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;               % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;               % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;               % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;           % do we vary with critical task saturation?
                Sim.doHF_v0     = true;               % do we have HF affect v0
                Sim.doHF_T      = true;               % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);
             
		   case 'scenBaseR_MVD_TD_PESV_BAST_H_I=10.0'
                Sim.mu = 0.40;
                Sim.ri = 10.0;
                Sim.doHF_TD     = true; 
                Sim.doHF_perc_s = true;      % do we have HF affect perception of distance gaps
                Sim.doHF_perc_v = true;      % do we have HF affect perception of speed differences
                Sim.doHF_tau    = true;      % do we have HF affect reaction time (attention lag)
                Sim.doHeterogeneity = true;  % do we vary with critical task saturation?
                Sim.doHF_v0     = true;      % do we have HF affect v0
                Sim.doHF_T      = true;     % do we have HF affect T
				inParPlot(Sim.doHF, Sim.doDisturbance, Sim.doHeterogeneity, Sim.doHF_tau, Sim.doHF_perc_s, Sim.doHF_perc_v,Sim.doHF_v0, Sim.doHF_T, Sim.doHF_TD);


        end
   
     % close current figs
        close all;
        
        % run simulations--------------------------------------------------
        Sim = simulation_fcn(Sim);
        %------------------------------------------------------------------
        % recieve information from the detectors---------------------------
        Sim = vLoopDet(Sim);
        %------------------------------------------------------------------
        % Calculate KPI's--------------------------------------------------
        % % My code goes here. Compute the followings:
        % % st     : stands for space - time
        % % TDC    : Total Distance Covered
        % % TTS    : Total Time Spent
        % % Edie.q : q = d(wmega) / |wmega|
        % % Edie.k : k = t(wmega) / |wmega|
        Sim = Edie(Sim);
        %% Display here the KPI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n---------KPI---------\n');
%         fprintf('TTS = %.3f [min]\n',Sim.Edie.TTS * 60);
%         fprintf('TTC_CRIT = %.1f  [s]\n',Sim.TTC_crit);
%         fprintf('TDC = %.1f  [km]\n',Sim.Edie.TDC);
%         fprintf('MEAN_SPEED = %.1f  [km/h]\n', Sim.Edie.u);
%         fprintf('CAPACITY = %.1f [veh/h]\n', Sim.Edie.q);
%         fprintf('DENSITY = %.1f [veh/km]\n', Sim.Edie.k)
%         fprintf('-------------------------\n');
        
        disp(table([Sim.Edie.TTS * 60; Sim.TTC_crit; Sim.Edie.u; Sim.Edie.q...
            ; Sim.Edie.k ],'VariableNames',{'Value'}, 'RowName',{'TTS', 'TTC_CRIT',...
            'MEAN_SPEED','CAPACITY', 'DENSITY' }));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figures----------------------------------------------------------
        if dofig
            % visualisations
            figs = visualisation_fcn(Sim);
        end
        %------------------------------------------------------------------
        % save data?-------------------------------------------------------
        if dosavedata
            disp('saving everything ...');
            
            if ~isdir(fullfile(pwd,'RainIntResults'))
                mkdir(pwd,'RainIntResults')
            end;
            fname = sprintf('RainIntResults/%s_%s.mat',Sim.scenario, tau_str);
            save(fname,'Sim');
        end;
        %------------------------------------------------------------------
        % save figures?----------------------------------------------------
        if dosavefig
            % we need the figs ...
            if ~dofig
                figs = visualisation_fcn(Sim);
            end;
            % export graphs to png
            if ~isdir(fullfile(pwd,'RainIntImages'))
                mkdir(pwd,'RainIntImages');
            end;
            exportgraphics(figs.Traj  ,sprintf('RainIntImages/%s_%s_traj.png',Sim.scenario,tau_str),'Resolution', 300);
            exportgraphics(figs.HF    ,sprintf('RainIntImages/%s_%s_hf.png',Sim.scenario,tau_str),'Resolution', 300);
            exportgraphics(figs.FDTD  ,sprintf('RainIntImages/%s_%s_fdtd.png',Sim.scenario,tau_str),'Resolution', 300);
            exportgraphics(figs.FD    ,sprintf('RainIntImages/%s_%s_fd.png',Sim.scenario,tau_str),'Resolution', 300);
            exportgraphics(figs.FDsel ,sprintf('RainIntImages/%s_%s_fdsel.png',Sim.scenario,tau_str),'Resolution', 300);
            exportgraphics(figs.Demand,sprintf('RainIntImages/%s_%s_demand.png',Sim.scenario,tau_str),'Resolution', 300);
			exportgraphics(figs.HEADACC,sprintf('RainIntImages/%s_%s_accheadspeed.png',Sim.scenario,tau_str),'Resolution', 300);
% 			exportgraphics(figs.PPD,sprintf('RainIntImages/%s_%s_virtualDetectors.png',Sim.scenario,tau_str),'Resolution', 300);
        end
        %------------------------------------------------------------------
        %% Here we clean the memory from the variable Sim, so as to 
        ... safeguard that every simulation contains unique parameters!
        clear Sim.Qvlp;
        clear Sim.Vvlp;
        clear Sim.Rhovlp;
        
        disp('DONE');
    end
end

Sim.numscen = numscen;

return;
%% Below follow cells with code to generate the two-column result graphs 
%  and the results table used in the manuscript

%% Compute TTS and num inc table
tbl_TTS = zeros(5, numscen);
tbl_INC = zeros(5, numscen);
tbl_TDC = zeros(5, numscen); % density
taus = [1,2,3,4];
for itau = 1:numel(taus)
    aTaustr = sprintf('tau=%.2d',taus(itau));
    for iscenario = 1:numscen
        aScenario = scenarios{iscenario};
        fname = sprintf('RainResults/%s_%s.mat',aScenario, aTaustr);
        load(fname);
        tbl_TTS(itau, iscenario) = Sim.TTS;
        tbl_INC(itau, iscenario) = Sim.numCollisions;
        tbl_TDC(itau, iscenario) = Sim.TDC;
    end
end;

fname = 'tbls2.xls'
xlswrite(fname,scenarios,'A1');
xlswrite(fname,[tbl_TTS; tbl_INC; tbl_TDC],'A2');

    
